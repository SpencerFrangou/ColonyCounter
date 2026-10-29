"""Run the colony formation analysis from a filled plate map.

    python cfa_run.py PlateMap.xlsx --out Results

Reads the plate map (see cfa_platemap.py), validates it, analyses every labelled well, and for each
comparison chooses the statistical model from the layout it finds:
  * blocked        treatment and control share plates          -> plate entered as a fixed blocking factor
  * partly blocked only some plates carry both                 -> plate factor, note on partial blocking
  * confounded     no plate carries both (one condition/plate) -> no plate term; confounding flagged
  * n < 3 wells in either group                                -> estimates only, no p-values
Count metrics: quasi-Poisson GLM with log(analysed area) offset -> rate ratio. Ratio metrics: OLS on log -> fold.
Linear metrics: OLS -> difference. Primary endpoints (colony number, stained area) Holm-adjusted per line;
secondary endpoints BH-adjusted per line. Games-Howell / Welch t, Mann-Whitney, Hedges g also reported.
Outputs: per_well_metrics.csv, colonies.csv, statistics_contrasts.csv, statistics_design.csv, Figure_<line>.(png|svg|pdf),
QC overlays with condition labels, and report.html.
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, os, sys, base64, itertools, datetime
import numpy as np, pandas as pd, cv2
from scipy import stats
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfa_pipeline as cp
from cfa_platemap import read_platemap, validate
from cfa_analysis import set_pub_style, OKABE_ITO, MM, hedges_g, stars

METRICS = {"whole_well_equiv": ("Colony number", "primary", "count", "Colonies per well\n(whole-well equiv.)"),
           "area_fraction": ("Stained area fraction", "primary", "log", "Stained area\nfraction"),
           "integrated_od_per_mm2": ("Integrated OD", "secondary", "log", "Integrated OD\n(per mm²)"),
           "mean_diam_mm": ("Mean colony diameter", "secondary", "linear", "Mean colony\ndiameter (mm)"),
           "mean_colony_area_mm2": ("Mean colony area", "secondary", "log", "Mean colony\narea (mm²)"),
           "mean_colony_od": ("Mean colony OD", "secondary", "linear", "Mean colony OD")}

# ----------------------------------------------------------------------------- image analysis
def analyse_all(W, settings, out):
    cp.RIM = float(settings.get("rim_fraction", 0.88)); cp.MIN_DIAM_MM = float(settings.get("min_diam_mm", 0.15))
    cp.GROWTH_AREA_MM2 = float(settings.get("growth_area_mm2", 950)); T = float(settings.get("threshold_od", 0.10))
    Cs, Ws, geo = [], [], []
    os.makedirs(f"{out}/QC_overlays", exist_ok=True)
    for plate, g in W.groupby("plate", sort=False):
        path = g.path.iloc[0]; im = cv2.imread(path)
        if im is None: sys.exit(f"cannot read {path}")
        wells, info = cp.detect_wells(im, dpi=cp.image_dpi(path)); byname = {w["name"]: w for w in wells}
        ov = im.copy()
        for r in g.itertuples():
            if r.well not in byname: print(f"warning: well {r.well} not detected on {plate}"); continue
            w = byname[r.well]
            C, wd = cp.analyse_well(im, w, info["px_per_mm"], T=T, overlay=ov, plate=plate)
            C["line"] = r.line; C["condition"] = r.condition; C["well"] = r.well
            wd.update(line=r.line, condition=r.condition, plate=plate, scan_date=r.scan_date, cells_seeded=r.cells_seeded, dpi=info["dpi"], px_per_mm_pitch=info["px_per_mm_pitch"])
            if r.cells_seeded: wd["colonies_per_100_seeded"] = wd["whole_well_equiv"] / float(r.cells_seeded) * 100
            cv2.putText(ov, f"{r.line}: {r.condition}", (int(w["cx"]) - 90, int(w["cy"] + w["r"] * 0.75)), cv2.FONT_HERSHEY_SIMPLEX, 0.9, (0, 0, 200), 2)
            Cs.append(C); Ws.append(wd)
        for w in wells:
            if w["name"] not in set(g.well): cv2.putText(ov, "not used", (int(w["cx"]) - 60, int(w["cy"])), cv2.FONT_HERSHEY_SIMPLEX, 0.9, (120, 120, 120), 2)
        xs = [w["cx"] for w in wells]; ys = [w["cy"] for w in wells]; R = info["r_ring_px"]
        y0, y1 = max(int(min(ys) - 1.3 * R), 0), int(max(ys) + 1.3 * R); x0, x1 = max(int(min(xs) - 1.3 * R), 0), int(max(xs) + 1.3 * R)
        cv2.imwrite(f"{out}/QC_overlays/{plate}_overlay.jpg", ov[y0:y1, x0:x1], [cv2.IMWRITE_JPEG_QUALITY, 85])
        geo.append(pd.DataFrame(wells).assign(plate=plate, px_per_mm=info["px_per_mm"]))
        print(f"{plate}: {len(g)} wells analysed, {info['px_per_mm']:.2f} px/mm")
    C = pd.concat(Cs, ignore_index=True); P = pd.DataFrame(Ws); G = pd.concat(geo, ignore_index=True)
    C.to_csv(f"{out}/colonies.csv", index=False); P.to_csv(f"{out}/per_well_metrics.csv", index=False); G.to_csv(f"{out}/well_geometry.csv", index=False)
    return C, P

# ----------------------------------------------------------------------------- design-aware statistics
def comparisons_for(W, controls, comps, settings):
    out = []
    if not comps.empty:
        for r in comps.itertuples(): out.append(dict(name=r.name or f"{r.treatment} vs {r.control}", line=r.line, treatment=r.treatment, control=r.control, source="Comparisons sheet"))
    lines_done = set(c["line"] for c in out)
    for line, g in W.groupby("line"):
        if line in lines_done: continue
        ctl = controls.get(line)
        for c in sorted(set(g.condition)):
            if c != ctl: out.append(dict(name=f"{c} vs {ctl}", line=line, treatment=c, control=ctl, source="Controls sheet"))
    if str(settings.get("all_pairs", "no")).lower().startswith("y"):
        for line, g in W.groupby("line"):
            have = {(c["treatment"], c["control"]) for c in out} | {(c["control"], c["treatment"]) for c in out}
            for a, b in itertools.combinations(sorted(set(g.condition)), 2):
                if (a, b) not in have: out.append(dict(name=f"{b} vs {a}", line=line, treatment=b, control=a, source="all pairs"))
    return out

def classify_design(P, line, t, c):
    pt = set(P[(P.line == line) & (P.condition == t)].plate); pc = set(P[(P.line == line) & (P.condition == c)].plate)
    both = pt & pc
    if not both: return "confounded", "no plate carries both conditions: plate and condition cannot be separated; plate-level staining/scanning differences are included in the estimate"
    if both == pt == pc: return "blocked", "every plate carries both conditions: plate entered as a blocking factor, which removes plate-level staining and scanning differences"
    return "partly blocked", f"{len(both)} plate(s) carry both conditions, others carry one: plate entered as a blocking factor; the estimate uses within-plate information where available"

def fit_contrast(d, col, kind, blocked):
    f = "y ~ cond" + (" + C(plate)" if blocked else "")
    if kind == "count":
        d = d.assign(y=d.n_colonies.astype(int), off=np.log(d.analysed_area_mm2))
        m = smf.glm(f, d, family=sm.families.Poisson(), offset=d.off).fit(scale="X2", use_t=True); model = f"quasi-Poisson GLM (phi={m.scale:.2f}), offset log(area)"
    elif kind == "log":
        d = d.assign(y=np.log(d[col].astype(float))); m = smf.ols(f, d).fit(); model = "OLS on log scale"
    else:
        d = d.assign(y=d[col].astype(float)); m = smf.ols(f, d).fit(); model = "OLS"
    if blocked: model += " + plate block"
    name = [n for n in m.params.index if n.startswith("cond[T.")][0]
    beta = m.params[name]; ci = m.conf_int().loc[name]; p = m.pvalues[name]
    if kind in ("count", "log"): return np.exp(beta), np.exp(ci[0]), np.exp(ci[1]), p, model, m
    return beta, ci[0], ci[1], p, model, m

def run_stats(P, comps, out):
    rows, design = [], []
    for cmp in comps:
        line, t, c = cmp["line"], cmp["treatment"], cmp["control"]
        d = P[(P.line == line) & (P.condition.isin([t, c]))].copy(); d["cond"] = pd.Categorical(d.condition, categories=[c, t])
        nt, nc = (d.condition == t).sum(), (d.condition == c).sum()
        dclass, dnote = classify_design(P, line, t, c); blocked = dclass != "confounded"
        testable = min(nt, nc) >= 3
        design.append(dict(comparison=cmp["name"], line=line, treatment=t, control=c, n_treatment=nt, n_control=nc,
                           plates_treatment=";".join(sorted(set(d[d.condition == t].plate))), plates_control=";".join(sorted(set(d[d.condition == c].plate))),
                           design=dclass, plate_term=("yes" if blocked else "no"), tests=("yes" if testable else "no: fewer than 3 wells in a group, estimates only"), note=dnote, source=cmp["source"]))
        for col, (label, family, kind, _) in METRICS.items():
            a = d[d.condition == c][col].astype(float).values; b = d[d.condition == t][col].astype(float).values
            try: est, lo, hi, p, model, m = fit_contrast(d, col, kind, blocked)
            except Exception as e: est = lo = hi = p = np.nan; model = f"not fitted ({e})"; m = None
            resid_p = np.nan
            if m is not None and len(d) >= 5:
                try: resid_p = stats.shapiro(m.resid_response if hasattr(m, "resid_response") else m.resid).pvalue
                except Exception: pass
            g, glo, ghi = hedges_g(a, b) if min(nt, nc) >= 2 else (np.nan, np.nan, np.nan)
            welch = stats.ttest_ind(b, a, equal_var=False).pvalue if testable else np.nan
            mwu = stats.mannwhitneyu(b, a).pvalue if testable else np.nan
            rows.append(dict(comparison=cmp["name"], line=line, metric=label, family=family, treatment=t, control=c, n_treatment=nt, n_control=nc, design=dclass,
                             model=model, estimate_type=("rate ratio" if kind == "count" else "fold change" if kind == "log" else "difference"),
                             mean_control=a.mean(), sd_control=a.std(ddof=1) if nc > 1 else np.nan, mean_treatment=b.mean(), sd_treatment=b.std(ddof=1) if nt > 1 else np.nan,
                             estimate=est, ci_low=lo, ci_high=hi, p_model=(p if testable else np.nan), p_welch=welch, p_mannwhitney=mwu,
                             shapiro_resid_p=resid_p, rank_test_recommended=("yes" if (resid_p == resid_p and resid_p < 0.01) else "no"),
                             hedges_g=g, g_ci_low=glo, g_ci_high=ghi))
    R = pd.DataFrame(rows); R["p_adj"] = np.nan; R["adjustment"] = ""
    for line in R.line.unique():
        for fam, meth, lab in [("primary", "holm", "Holm"), ("secondary", "fdr_bh", "BH q")]:
            idx = (R.line == line) & (R.family == fam) & R.p_model.notna()
            if idx.any(): R.loc[idx, "p_adj"] = multipletests(R.loc[idx, "p_model"], method=meth)[1]; R.loc[idx, "adjustment"] = f"{lab} ({int(idx.sum())} tests, line {line})"
    D = pd.DataFrame(design)
    R.to_csv(f"{out}/statistics_contrasts.csv", index=False); D.to_csv(f"{out}/statistics_design.csv", index=False)
    return R, D

# ----------------------------------------------------------------------------- figures
def figure_line(P, C, line, R, comps, out):
    set_pub_style(); plt.rcParams["axes.unicode_minus"] = False
    conds = list(dict.fromkeys(list(P[P.line == line].condition)))
    ctl_of = {c["treatment"]: c["control"] for c in comps if c["line"] == line}
    cols = {c: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(conds)}; k = len(conds)
    panels = ["whole_well_equiv", "area_fraction", "integrated_od_per_mm2", "mean_diam_mm"]
    ncol = 2 if k <= 5 else 1; nrow = int(np.ceil(len(panels) / ncol)); ph = 55 * MM
    Rl = R[R.line == line]; nF = len(Rl[Rl.metric.isin(["Colony number", "Stained area fraction", "Integrated OD", "Mean colony area"])])
    bottom_h = ph if ncol == 2 else ph + max(45, 4.5 * nF + 12) * MM
    fig = plt.figure(figsize=(180 * MM, nrow * ph + bottom_h + 16 * MM)); top_frac = nrow * ph / (nrow * ph + bottom_h + 16 * MM)
    gs = fig.add_gridspec(nrow, ncol, left=0.11, right=0.98, top=0.96, bottom=0.96 - top_frac + 0.06, hspace=0.6, wspace=0.5)
    if ncol == 2:
        gsb = fig.add_gridspec(1, 2, left=0.11, right=0.98, top=0.96 - top_frac - 0.02, bottom=0.08, wspace=1.3, width_ratios=[1, 0.9]); axE, axF = fig.add_subplot(gsb[0, 0]), fig.add_subplot(gsb[0, 1])
    else:
        gsb = fig.add_gridspec(2, 1, left=0.32, right=0.98, top=0.96 - top_frac - 0.02, bottom=0.05, hspace=0.5, height_ratios=[ph, bottom_h - ph]); axE, axF = fig.add_subplot(gsb[0, 0]), fig.add_subplot(gsb[1, 0])
    rng = np.random.default_rng(0); labx = -0.26 if ncol == 2 else -0.1
    for pi, col in enumerate(panels):
        ax = fig.add_subplot(gs[pi // ncol, pi % ncol]); label, family, kind, ylab = METRICS[col]
        groups = [P[(P.line == line) & (P.condition == c)][col].astype(float).values for c in conds]
        for i, (c, v) in enumerate(zip(conds, groups)):
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=cols[c], edgecolor="white", linewidth=0.4, zorder=3)
            if len(v) > 1: ax.errorbar(i + 0.32, v.mean(), yerr=v.std(ddof=1), fmt="_", color="#222", capsize=2.5, markersize=8, elinewidth=0.9, capthick=0.9, zorder=4)
        xl = [c.replace("-", "-\n") if len(c) > 8 else c for c in conds]; rot = 40 if (k > 4 and max(len(x) for x in xl) > 6) else 0
        ax.set_xticks(range(k)); ax.set_xticklabels(xl, rotation=rot, ha=("right" if rot else "center")); ax.set_xlim(-0.6, k - 0.2); ax.set_ylabel(ylab)
        ymax = max(v.max() for v in groups); ymin = min(v.min() for v in groups); yr = max(ymax - ymin, 1e-9)
        ax.set_ylim(0 if col == "whole_well_equiv" else ymin - 0.1 * yr, ymax + 0.32 * yr)
        for j, c in enumerate(conds):
            if c not in ctl_of: continue
            rr = Rl[(Rl.metric == label) & (Rl.treatment == c) & (Rl.control == ctl_of[c])]
            if len(rr) and rr.p_adj.notna().iloc[0]: ax.text(j, ymax + 0.08 * yr, stars(rr.p_adj.iloc[0]), ha="center", va="bottom", fontsize=8)
        ax.text(labx, 1.06, "ABCD"[pi], transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    Cin = C[(C.fully_inside) & (C.line == line)]
    for c in conds:
        v = np.sort(Cin[Cin.condition == c].eq_diam_mm.values)
        if len(v): axE.step(v, np.arange(1, len(v) + 1) / len(v), color=cols[c], lw=1.0, label=f"{c} (n = {len(v)})")
    axE.set_xlabel("Colony equivalent\ndiameter (mm)"); axE.set_ylabel("Cumulative\nfraction"); axE.legend(frameon=False, loc="lower right", fontsize=8, handlelength=1.2)
    axE.text(labx if ncol == 2 else -0.14, 1.06, "E", transform=axE.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    rr = Rl[Rl.metric.isin(["Colony number", "Stained area fraction", "Integrated OD", "Mean colony area"]) & Rl.estimate.notna()].copy()
    ms = {"Colony number": "Colony no.", "Stained area fraction": "Stained area", "Integrated OD": "Integrated OD", "Mean colony area": "Colony area"}
    rr["lab"] = rr.metric.map(ms) + ", " + rr.treatment + " vs " + rr.control; rr = rr.iloc[::-1]; y = np.arange(len(rr))
    for yi, (_, q) in zip(y, rr.iterrows()):
        axF.plot([q.ci_low, q.ci_high], [yi, yi], color=cols.get(q.treatment, "#333"), lw=1.0); axF.plot(q.estimate, yi, "o", color=cols.get(q.treatment, "#333"), ms=4)
    axF.axvline(1, color="#888", lw=0.75, ls="--"); axF.set_xscale("log"); axF.set_yticks(y); axF.set_yticklabels(rr.lab, fontsize=8)
    if len(rr):
        lo_, hi_ = rr.ci_low.min(), rr.ci_high.max(); ticks = [t for t in [0.2, 0.25, 0.33, 0.5, 0.67, 1, 1.5, 2, 3, 4, 5] if lo_ * 0.9 <= t <= hi_ * 1.1]
        if len(ticks) < 3: ticks = sorted(set(np.round(np.exp(np.linspace(np.log(lo_ * 0.95), np.log(hi_ * 1.05), 5)), 2)) | {1.0})
        axF.xaxis.set_major_locator(FixedLocator(ticks)); axF.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}")); axF.xaxis.set_minor_locator(NullLocator())
    axF.set_xlabel("Ratio, treatment / control\n(95% CI)"); axF.set_ylim(-0.7, max(len(rr) - 0.3, 0.3))
    axF.text(-0.8 if ncol == 2 else -0.3, 1.04, "F", transform=axF.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    with matplotlib.rc_context({"pdf.use14corefonts": True}): fig.savefig(f"{out}/Figure_{line}.pdf")
    fig.savefig(f"{out}/Figure_{line}.svg"); fig.savefig(f"{out}/Figure_{line}.png", dpi=600); plt.close(fig)

# ----------------------------------------------------------------------------- report
def report(P, R, D, settings, out, platemap):
    def img(p): return "data:image/jpeg;base64," + base64.b64encode(open(p, "rb").read()).decode()
    summ = P.groupby(["line", "condition"]).agg(wells=("well", "size"), plates=("plate", "nunique"), colonies_per_well=("whole_well_equiv", "mean"), sd=("whole_well_equiv", "std"),
                                                 area_fraction=("area_fraction", "mean"), integrated_od_per_mm2=("integrated_od_per_mm2", "mean"), mean_diam_mm=("mean_diam_mm", "mean")).round(3).reset_index()
    cols = ["comparison", "line", "metric", "design", "model", "estimate_type", "estimate", "ci_low", "ci_high", "p_model", "p_adj", "adjustment", "p_welch", "p_mannwhitney", "rank_test_recommended", "hedges_g"]
    h = [f"<html><head><meta charset='utf-8'><title>Colony formation analysis report</title><style>body{{font-family:Helvetica,Arial,sans-serif;font-size:13px;margin:30px}}table{{border-collapse:collapse}}td,th{{border:1px solid #ccc;padding:3px 6px;text-align:right}}th{{background:#eee}}td:first-child,th:first-child{{text-align:left}}img{{max-width:900px;display:block;margin:6px 0 18px}}h2{{margin-top:30px}}</style></head><body>",
         f"<h1>Colony formation analysis report</h1><p>Plate map: {os.path.basename(platemap)} &middot; run {datetime.datetime.now():%Y-%m-%d %H:%M} &middot; settings: " + ", ".join(f"{k} = {v}" for k, v in settings.items()) + "</p>",
         "<h2>1. Design decisions</h2><p>How each comparison was analysed, decided from the plate layout.</p>", D.to_html(index=False, border=0),
         "<h2>2. Summary per condition</h2>", summ.to_html(index=False, border=0),
         "<h2>3. Statistics</h2><p>Primary endpoints (colony number, stained area fraction) are Holm-adjusted within each cell line; secondary endpoints carry Benjamini-Hochberg q-values. Where the residual Shapiro-Wilk p is below 0.01 the rank-based Mann-Whitney p is the safer one to quote.</p>",
         R[cols].round(4).to_html(index=False, border=0)]
    for line in P.line.unique():
        h.append(f"<h2>Figure: {line}</h2><img src='{img(f'{out}/Figure_{line}.png')}'>")
    h.append("<h2>QC overlays</h2><p>Blue = fitted rim; red = analysed disk; green = counted colonies; label = cell line and condition from the plate map. Check every well.</p>")
    for p in sorted(os.listdir(f"{out}/QC_overlays")): h.append(f"<h3>{p}</h3><img src='{img(f'{out}/QC_overlays/{p}')}'>")
    h.append("</body></html>"); open(f"{out}/report.html", "w").write("\n".join(h))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("platemap"); ap.add_argument("--out", default="Results"); a = ap.parse_args()
    W, controls, comps_df, settings = read_platemap(a.platemap); problems = validate(W, controls, comps_df)
    if problems: sys.exit("Plate map problems:\n  " + "\n  ".join(problems))
    os.makedirs(a.out, exist_ok=True); print(f"{len(W)} wells on {W.plate.nunique()} plates; lines: {', '.join(W.line.unique())}")
    C, P = analyse_all(W, settings, a.out)
    comps = comparisons_for(W, controls, comps_df, settings); print(f"{len(comps)} comparisons")
    R, D = run_stats(P, comps, a.out)
    for line in P.line.unique(): figure_line(P, C, line, R, comps, a.out)
    report(P, R, D, settings, a.out, a.platemap)
    pd.set_option("display.width", 250); print(D[["comparison", "n_treatment", "n_control", "design", "plate_term", "tests"]].to_string(index=False))
    print(R[R.family == "primary"][["comparison", "metric", "estimate", "ci_low", "ci_high", "p_adj", "hedges_g"]].round(3).to_string(index=False))
    print(f"report: {a.out}/report.html")

if __name__ == "__main__": main()
