"""Colony formation assay: batch quantification, statistics and publication figures.

Usage:
    python cfa_analysis.py samples.csv --out results/ [--threshold 0.10] [--all-pairs]

samples.csv columns:
    file        path to the scanned plate image (jpg/png/tif)
    line        cell line / experiment group; statistics and figures are produced per line
    condition   condition label (one plate = one condition; six wells = replicates)
    reference   1 for the control condition of that line (exactly one per line), else 0
    order       (optional) integer display order within the line
    cells_seeded (optional) cells seeded per well; if given, colonies per 100 seeded cells is reported

Design handled: any number of conditions per line (2 ... n). Contrasts are condition-vs-reference
by default (many-to-one); --all-pairs tests every pair. Primary endpoints (colony number, stained
area fraction): Holm within the primary family per line. Secondary endpoints: Benjamini-Hochberg
q-values per line. Count metrics: quasi-Poisson GLM with log(analysed area) offset -> rate ratios.
Ratio metrics: OLS on log scale -> fold changes. Linear metrics: OLS -> differences.
Post hoc: Games-Howell; effect size: Hedges' g with bootstrap CI.
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, itertools, os, sys
import numpy as np, pandas as pd, cv2
from scipy import stats
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfa_pipeline as cp

# ----------------------------------------------------------------------------- figure style
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#000000", "#999999", "#882255"]
MM = 1 / 25.4
def set_pub_style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 10, "axes.titlesize": 10, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.75, "lines.linewidth": 1.0, "xtick.major.width": 0.75, "ytick.major.width": 0.75,
        "xtick.major.size": 3, "ytick.major.size": 3, "axes.spines.top": False, "axes.spines.right": False,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 100})

# ----------------------------------------------------------------------------- statistics helpers
def hedges_g(a, b, n_boot=4000, seed=0):
    na, nb = len(a), len(b); J = 1 - 3 / (4 * (na + nb) - 9)
    def g_of(x, y):
        s = np.sqrt(((len(x) - 1) * x.var(ddof=1) + (len(y) - 1) * y.var(ddof=1)) / (len(x) + len(y) - 2))
        return J * (y.mean() - x.mean()) / s if s > 0 else np.nan
    rng = np.random.default_rng(seed)
    gs = [g_of(rng.choice(a, na), rng.choice(b, nb)) for _ in range(n_boot)]
    return g_of(a, b), np.nanpercentile(gs, 2.5), np.nanpercentile(gs, 97.5)

def games_howell_p(a, b, k):
    na, nb = len(a), len(b); va, vb = a.var(ddof=1) / na, b.var(ddof=1) / nb
    t = (b.mean() - a.mean()) / np.sqrt(va + vb); df = (va + vb) ** 2 / (va ** 2 / (na - 1) + vb ** 2 / (nb - 1))
    return float(stats.studentized_range.sf(abs(t) * np.sqrt(2), k, df))

METRICS = {  # column: (label, family, model kind, axis label)
    "whole_well_equiv": ("Colony number", "primary", "count", "Colonies per well\n(whole-well equivalent)"),
    "area_fraction": ("Stained area fraction", "primary", "log", "Stained area fraction"),
    "integrated_od_per_mm2": ("Integrated OD", "secondary", "log", "Integrated OD (mm$^{-2}$)"),
    "mean_diam_mm": ("Mean colony diameter", "secondary", "linear", "Mean colony diameter (mm)"),
    "mean_colony_area_mm2": ("Mean colony area", "secondary", "log", "Mean colony area (mm$^2$)"),
    "mean_colony_od": ("Mean colony OD", "secondary", "linear", "Mean colony OD"),
}

def run_stats(P, line, conds, ref, all_pairs=False):
    sub = P[P.plate.isin(conds)].copy(); sub["cond"] = pd.Categorical(sub.plate, categories=[ref] + [c for c in conds if c != ref])
    k = len(conds); rows, omni = [], []
    pairs = list(itertools.combinations(conds, 2)) if all_pairs else [(ref, c) for c in conds if c != ref]
    for col, (label, family, kind, _) in METRICS.items():
        groups = {c: sub[sub.plate == c][col].astype(float).values for c in conds}
        if k >= 3:
            W = stats.alexandergovern(*groups.values()).pvalue; K = stats.kruskal(*groups.values()).pvalue
        else:
            W = stats.ttest_ind(groups[conds[0]], groups[conds[1]], equal_var=False).pvalue; K = stats.mannwhitneyu(groups[conds[0]], groups[conds[1]]).pvalue
        Lv = stats.levene(*groups.values()).pvalue; Sh = stats.shapiro(np.concatenate([g - g.mean() for g in groups.values()])).pvalue
        omni.append(dict(line=line, metric=label, family=family, omnibus=("Welch ANOVA" if k >= 3 else "Welch t"), omnibus_p=W,
                         rank_test=("Kruskal-Wallis" if k >= 3 else "Mann-Whitney"), rank_p=K, levene_p=Lv, shapiro_resid_p=Sh))
        d = sub.copy()
        if kind == "count":
            d["y"] = d["n_colonies"].astype(int); d["off"] = np.log(d["analysed_area_mm2"])
            m = smf.glm("y ~ cond", d, family=sm.families.Poisson(), offset=d["off"]).fit(scale="X2", use_t=True); model = f"quasi-Poisson GLM (phi = {m.scale:.2f}), offset log(area)"
        elif kind == "log":
            d["y"] = np.log(d[col].astype(float)); m = smf.ols("y ~ cond", d).fit(); model = "OLS on log scale"
        else:
            d["y"] = d[col].astype(float); m = smf.ols("y ~ cond", d).fit(); model = "OLS"
        for a, b in pairs:
            L = np.array([(1.0 if n == f"cond[T.{b}]" else (-1.0 if n == f"cond[T.{a}]" else 0.0)) for n in m.params.index])
            tt = m.t_test(L); beta = float(np.ravel(tt.effect)[0]); ci = np.ravel(tt.conf_int()); p = float(np.ravel(tt.pvalue)[0])
            if kind in ("count", "log"): est, lo, hi, et = np.exp(beta), np.exp(ci[0]), np.exp(ci[1]), ("rate ratio" if kind == "count" else "fold change")
            else: est, lo, hi, et = beta, ci[0], ci[1], "difference"
            g, glo, ghi = hedges_g(groups[a], groups[b])
            rows.append(dict(line=line, metric=label, family=family, model=model, reference=a, condition=b, contrast=f"{b} vs {a}",
                             n_ref=len(groups[a]), n_cond=len(groups[b]), mean_ref=groups[a].mean(), sd_ref=groups[a].std(ddof=1),
                             mean_cond=groups[b].mean(), sd_cond=groups[b].std(ddof=1), estimate_type=et, estimate=est, ci_low=lo, ci_high=hi,
                             p_model=p, p_games_howell=games_howell_p(groups[a], groups[b], k), hedges_g=g, g_ci_low=glo, g_ci_high=ghi))
    R = pd.DataFrame(rows); R["p_adj"] = np.nan; R["adjustment"] = ""
    pr = R.family == "primary"; se = R.family == "secondary"
    R.loc[pr, "p_adj"] = multipletests(R.loc[pr, "p_model"], method="holm")[1]; R.loc[pr, "adjustment"] = f"Holm (primary family, {pr.sum()} tests)"
    R.loc[se, "p_adj"] = multipletests(R.loc[se, "p_model"], method="fdr_bh")[1]; R.loc[se, "adjustment"] = f"BH FDR q (secondary family, {se.sum()} tests)"
    R["q_bh_all"] = multipletests(R.p_model, method="fdr_bh")[1]
    return R, pd.DataFrame(omni)

# ----------------------------------------------------------------------------- figure
def stars(p): return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

def pub_figure(P, C, line, conds, ref, R, out_prefix, all_pairs=False):
    set_pub_style(); k = len(conds); short = {c: (c.split("-", 1)[1] if "-" in c else c) for c in conds}
    cols = {c: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(conds)}
    panels = ["whole_well_equiv", "area_fraction", "integrated_od_per_mm2", "mean_diam_mm"]
    ncol = 2 if k <= 5 else 1
    nrow_dot = int(np.ceil(len(panels) / ncol))
    n_contr = len(R[R.metric.isin(["Colony number", "Stained area fraction", "Integrated OD", "Mean colony area"])])
    fig_w = 180 * MM; ph = 52 * MM
    bottom_h = ph if ncol == 2 else ph + max(45, 4.5 * n_contr + 12) * MM
    fig_h = nrow_dot * ph + bottom_h + 16 * MM
    fig = plt.figure(figsize=(fig_w, fig_h))
    top_frac = (nrow_dot * ph) / fig_h
    gs = fig.add_gridspec(nrow_dot, ncol, left=0.10, right=0.98, top=0.96, bottom=0.96 - top_frac + 0.06, hspace=0.55, wspace=0.45)
    if ncol == 2:
        gsb = fig.add_gridspec(1, 2, left=0.10, right=0.98, top=0.96 - top_frac - 0.02, bottom=0.07, wspace=1.25, width_ratios=[1.0, 0.9])
        ax_ecdf, ax_forest = fig.add_subplot(gsb[0, 0]), fig.add_subplot(gsb[0, 1])
    else:
        gsb = fig.add_gridspec(2, 1, left=0.30, right=0.98, top=0.96 - top_frac - 0.02, bottom=0.05, hspace=0.5, height_ratios=[ph, bottom_h - ph])
        ax_ecdf, ax_forest = fig.add_subplot(gsb[0, 0]), fig.add_subplot(gsb[1, 0])
    lab_x = -0.22 if ncol == 2 else -0.09
    axes = [fig.add_subplot(gs[i // ncol, i % ncol]) for i in range(len(panels))]
    rng = np.random.default_rng(0)
    for pi, (ax, col) in enumerate(zip(axes, panels)):
        label, family, kind, ylab = METRICS[col]
        groups = [P[P.plate == c][col].astype(float).values for c in conds]
        for i, (c, v) in enumerate(zip(conds, groups)):
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=14, color=cols[c], edgecolor="white", linewidth=0.4, zorder=3)
            ax.errorbar(i + 0.32, v.mean(), yerr=v.std(ddof=1), fmt="_", color="#222", capsize=2.5, markersize=8, elinewidth=0.9, capthick=0.9, zorder=4)
        xl = [short[c].replace("-", "-\n") if len(short[c]) > 8 else short[c] for c in conds]
        rot = 40 if (k > 4 and max(len(x) for x in xl) > 6) else 0
        ax.set_xticks(range(k)); ax.set_xticklabels(xl, rotation=rot, ha=("right" if rot else "center"))
        ax.set_xlim(-0.5, k - 0.5 + 0.5); ax.set_ylabel(ylab)
        ymax = max(v.max() for v in groups); ymin = min(v.min() for v in groups); yr = ymax - ymin
        if col in ("whole_well_equiv",): ax.set_ylim(0, ymax + 0.45 * yr)
        else: ax.set_ylim(ymin - 0.1 * yr, ymax + 0.45 * yr)
        # significance vs reference (or brackets for all pairs when k == 2)
        r = R[R.metric == label]; adjlab = "Holm p" if family == "primary" else "q"
        ri = conds.index(ref); y0 = ymax + 0.08 * yr
        if not all_pairs or k == 2:
            for j, c in enumerate(conds):
                if c == ref: continue
                rr = r[r.condition == c].iloc[0]; txt = stars(rr.p_adj)
                ax.text(j, y0, txt, ha="center", va="bottom", fontsize=(8 if txt != "ns" else 7), color="#222")
        ax.text(0.0, 1.0, f"n = {len(groups[0])} wells per condition", transform=ax.transAxes, fontsize=7, color="#555", ha="left", va="bottom")
        ax.text(lab_x, 1.06, "ABCDEFGH"[pi], transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    # bottom row: colony size ECDF (left) and effect-size forest plot (right)
    ax = ax_ecdf; Cin = C[C.fully_inside]
    for c in conds:
        v = np.sort(Cin[Cin.plate == c].eq_diam_mm.values); ax.step(v, np.arange(1, len(v) + 1) / len(v), color=cols[c], lw=1.0, label=f"{short[c]} (n = {len(v)})")
    ax.set_xlabel("Colony equivalent diameter (mm)"); ax.set_ylabel("Cumulative fraction"); ax.legend(frameon=False, loc="lower right", handlelength=1.2, fontsize=7)
    ax.text(lab_x if ncol == 2 else -0.13, 1.06, "ABCDEFGH"[len(panels)], transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    ax.text(0.0, 1.0, "colonies pooled across wells; descriptive", transform=ax.transAxes, fontsize=7, color="#555", ha="left", va="bottom")
    if True:
        ax = ax_forest
        rr = R[R.metric.isin(["Colony number", "Stained area fraction", "Integrated OD", "Mean colony area"])].copy()
        mshort = {"Colony number": "Colony no.", "Stained area fraction": "Stained area", "Integrated OD": "Integrated OD", "Mean colony area": "Colony area"}
        rr["lab"] = rr.metric.map(mshort) + ", " + rr.condition.map(short)
        rr = rr.iloc[::-1]; y = np.arange(len(rr))
        for yi, (_, q) in zip(y, rr.iterrows()):
            ax.plot([q.ci_low, q.ci_high], [yi, yi], color=cols[q.condition], lw=1.0); ax.plot(q.estimate, yi, "o", color=cols[q.condition], ms=4)
        ax.axvline(1, color="#888", lw=0.75, ls="--"); ax.set_xscale("log"); ax.set_yticks(y); ax.set_yticklabels(rr.lab, fontsize=7)
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator
        lo_, hi_ = rr.ci_low.min(), rr.ci_high.max(); ticks = [t for t in [0.25, 0.33, 0.5, 0.67, 1, 1.5, 2, 3, 4] if lo_ * 0.9 <= t <= hi_ * 1.1]
        ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v:g}")); ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel(f"Ratio vs {short[ref]} (95% CI)"); ax.set_ylim(-0.7, len(rr) - 0.3)
        ax.text(-0.75 if ncol == 2 else -0.3, 1.02, "ABCDEFGH"[len(panels) + 1], transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    fig.savefig(out_prefix + ".png", dpi=600); fig.savefig(out_prefix + ".svg"); fig.savefig(out_prefix + ".pdf"); plt.close(fig)
    # legend text
    prim = R[R.family == "primary"]; sec = R[R.family == "secondary"]
    leg = [f"Figure. Colony formation of {line} cells ({', '.join(short[c] for c in conds)}; reference = {short[ref]}).",
           f"(A) Colony number per well (whole-well equivalent, colonies counted within the central 0.88 R disk scaled to 950 mm2). (B) Stained area fraction. (C) Integrated optical density per mm2. (D) Mean colony equivalent diameter (colonies fully inside the analysed disk). Points are individual wells (n = 6 per condition, one plate per condition); bars show mean ± SD. (E) Cumulative distribution of colony diameter, colonies pooled across wells (descriptive). (F) Effect estimates with 95% CI: rate ratios from a quasi-Poisson GLM with log(analysed area) offset for colony number, fold changes from linear models on the log scale for area and intensity metrics.",
           f"Significance markers in A-D refer to each condition versus {short[ref]}: primary endpoints (A, B) Holm-adjusted p within the primary family ({len(prim)} tests); secondary endpoints (C, D) Benjamini-Hochberg FDR q ({len(sec)} tests). *** p or q < 0.001; ** < 0.01; * < 0.05; ns, not significant. Games-Howell post hoc p-values and Hedges' g with bootstrap 95% CI are given in the statistics table.",
           "Images: crystal violet stained 6-well plates scanned at 300 dpi; OD = log10(I_background / I_green); fixed segmentation threshold 0.10 OD; well rims fitted per well; outer 12% of the well radius excluded."]
    open(out_prefix + "_legend.txt", "w").write("\n\n".join(leg))

# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("samples"); ap.add_argument("--out", default="cfa_results"); ap.add_argument("--threshold", type=float, default=0.10)
    ap.add_argument("--all-pairs", action="store_true"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); S = pd.read_csv(a.samples)
    if "order" not in S: S["order"] = range(len(S))
    S = S.sort_values(["line", "order"]); Cs, Ws, Ls = [], [], []
    for r in S.itertuples():
        name = f"{r.line}-{r.condition}"
        C, W, L, info = cp.analyse_plate(r.file, name, T=a.threshold, out_overlay=f"{a.out}/{name}_overlay.jpg")
        C["line"] = r.line; C["condition"] = r.condition; W["line"] = r.line; W["condition"] = r.condition
        if "cells_seeded" in S and pd.notna(r.cells_seeded): W["colonies_per_100_seeded"] = W.whole_well_equiv / r.cells_seeded * 100
        print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in info.items()}); Cs.append(C); Ws.append(W); Ls.append(L)
    C = pd.concat(Cs, ignore_index=True); P = pd.concat(Ws, ignore_index=True); L = pd.concat(Ls, ignore_index=True)
    C.to_csv(f"{a.out}/colonies.csv", index=False); P.to_csv(f"{a.out}/per_well_metrics.csv", index=False); L.to_csv(f"{a.out}/well_geometry.csv", index=False)
    allR, allO = [], []
    for line, g in S.groupby("line", sort=False):
        conds = [f"{line}-{c}" for c in g.condition]; refs = g[g.reference == 1].condition.tolist()
        if len(refs) != 1: raise SystemExit(f"line {line}: exactly one reference condition required")
        ref = f"{line}-{refs[0]}"
        R, O = run_stats(P, line, conds, ref, all_pairs=a.all_pairs); allR.append(R); allO.append(O)
        pub_figure(P, C, line, conds, ref, R, f"{a.out}/Figure_{line}", all_pairs=a.all_pairs)
    R = pd.concat(allR, ignore_index=True); O = pd.concat(allO, ignore_index=True)
    R.to_csv(f"{a.out}/statistics_contrasts.csv", index=False); O.to_csv(f"{a.out}/statistics_omnibus.csv", index=False)
    pd.set_option("display.width", 250)
    print(R[["line", "metric", "contrast", "estimate_type", "estimate", "ci_low", "ci_high", "p_model", "p_adj", "p_games_howell", "hedges_g"]].round(4).to_string(index=False))

if __name__ == "__main__": main()
