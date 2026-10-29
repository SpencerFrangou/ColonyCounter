"""Re-derive the two empirical parameters after a change of scanner, resolution, plate type or stain protocol.

    python cfa_calibrate.py PlateMap.xlsx --out Calibration

1. Threshold sweep: colonies per well at OD thresholds 0.06 ... 0.25 for every plate. Choose the threshold at the
   low end of the plateau (counts change < 5 % per 0.02 OD step); 0.10 OD was the value for the 2026-10-29 scans.
2. Edge profile: density of segmented objects per mm2 in radial bands of the well (r/R from 0.50 to 1.00) at that
   threshold, pooled over all wells. Choose the rim fraction as the last band before the density rises; 0.88 was the
   value for the Epson V500 scans. Writes threshold_sweep.csv, edge_profile.csv and calibration.png.
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, os, sys, numpy as np, pandas as pd, cv2
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cfa_pipeline as cp
from cfa_platemap import read_platemap

ap = argparse.ArgumentParser(); ap.add_argument("platemap"); ap.add_argument("--out", default="Calibration"); ap.add_argument("--threshold", type=float, default=None)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
W, controls, comps, settings = read_platemap(a.platemap)
T0 = a.threshold or float(settings.get("threshold_od", 0.10)); Ts = [0.06, 0.08, 0.10, 0.12, 0.14, 0.17, 0.20, 0.25]
sweep, prof = [], []
for plate, g in W.groupby("plate", sort=False):
    im = cv2.imread(g.path.iloc[0]); wells, info = cp.detect_wells(im, dpi=cp.image_dpi(g.path.iloc[0])); byname = {w["name"]: w for w in wells}
    for r in g.itertuples():
        w = byname.get(r.well)
        if w is None: continue
        for T in Ts:
            _, wd = cp.analyse_well(im, w, info["px_per_mm"], T=T); sweep.append(dict(plate=plate, well=r.well, condition=r.condition, threshold=T, n_colonies=wd["n_colonies"]))
        # edge profile: segment the whole well (rim fraction 1.0) and bin colony centroids by r/R
        saved = cp.RIM; cp.RIM = 1.0; cp.BAND = 2.0
        C, _ = cp.analyse_well(im, w, info["px_per_mm"], T=T0); cp.RIM = saved; cp.BAND = 0.80
        R = w["r"]; edges = np.arange(0.50, 1.001, 0.04)
        for lo, hi in zip(edges[:-1], edges[1:]):
            area = np.pi * (hi ** 2 - lo ** 2) * (R / info["px_per_mm"]) ** 2
            n = int(((C.r_norm >= lo) & (C.r_norm < hi)).sum()) if len(C) else 0
            prof.append(dict(plate=plate, well=r.well, r_lo=lo, r_hi=hi, n=n, area_mm2=area))
    print(plate, "done")
S = pd.DataFrame(sweep); P = pd.DataFrame(prof); S.to_csv(f"{a.out}/threshold_sweep.csv", index=False); P.to_csv(f"{a.out}/edge_profile.csv", index=False)
st = S.groupby(["plate", "threshold"]).n_colonies.mean().unstack(); pr = P.groupby("r_hi").agg(n=("n", "sum"), area=("area_mm2", "sum")); pr["density"] = pr.n / pr.area
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for p in st.index: ax[0].plot(st.columns, st.loc[p], marker="o", ms=3, label=p)
ax[0].axvline(T0, color="#888", ls="--"); ax[0].set_xlabel("OD threshold"); ax[0].set_ylabel("mean colonies per well"); ax[0].legend(fontsize=7, frameon=False); ax[0].set_title("Threshold sweep")
ax[1].step(pr.index, pr.density, where="pre"); ax[1].axvline(float(settings.get("rim_fraction", 0.88)), color="#888", ls="--"); ax[1].set_xlabel("r / R (band upper edge)"); ax[1].set_ylabel("objects per mm2"); ax[1].set_title(f"Edge profile at {T0} OD")
fig.tight_layout(); fig.savefig(f"{a.out}/calibration.png", dpi=150)
print(st.round(1).to_string()); print(pr.round(3).to_string()); print("written to", a.out)
