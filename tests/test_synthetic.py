"""End-to-end check on synthetic plates: generated colony counts must be recovered within tolerance.
Run:  python tests/test_synthetic.py        (takes about a minute; needs no real data)"""
import os, subprocess, sys, tempfile, numpy as np, pandas as pd
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = tempfile.mkdtemp(prefix="colonycounter_test_")
subprocess.check_call([sys.executable, f"{root}/examples/make_synthetic_plates.py", "--out", out, "--seed", "7"])
subprocess.check_call([sys.executable, f"{root}/colonycounter/cfa_run.py", f"{out}/PlateMap.xlsx", "--out", f"{out}/Results"])
truth = {l.split("\t")[0]: list(map(int, l.strip().split("\t")[1:])) for l in open(f"{out}/ground_truth.txt")}
P = pd.read_csv(f"{out}/Results/per_well_metrics.csv")
ok = True
for name, counts in truth.items():
    got = P[P.plate == f"synthetic_{name}"].sort_values("well").whole_well_equiv.values
    err = np.abs(got - np.array(counts)) / np.array(counts)
    print(f"{name}: drawn {counts} recovered {np.round(got).astype(int).tolist()} max rel. error {err.max():.2f}")
    ok &= err.max() < 0.25
R = pd.read_csv(f"{out}/Results/statistics_contrasts.csv"); rr = R[R.metric == "Colony number"].estimate.iloc[0]
print(f"rate ratio treated/control {rr:.2f} (expected about 0.5)"); ok &= 0.35 < rr < 0.65
print("PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
