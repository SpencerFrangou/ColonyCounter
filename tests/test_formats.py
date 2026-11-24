"""Synthetic 12- and 24-well plates: every well is detected and the plate-level colony totals are recovered."""
import os, subprocess, sys, tempfile
import numpy as np, pandas as pd, pytest
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

@pytest.mark.parametrize("fmt", [12, 24])
def test_other_formats_are_detected_and_counted(fmt):
    out = tempfile.mkdtemp(prefix=f"colonycounter_{fmt}_")
    subprocess.check_call([sys.executable, f"{ROOT}/examples/make_synthetic_plates.py", "--out", out, "--seed", "3", "--format", str(fmt)])
    subprocess.check_call([sys.executable, f"{ROOT}/colonycounter/cfa_run.py", f"{out}/PlateMap.xlsx", "--out", f"{out}/Results"])
    truth = {l.split("\t")[0]: list(map(int, l.strip().split("\t")[1:])) for l in open(f"{out}/ground_truth.txt")}
    P = pd.read_csv(f"{out}/Results/per_well_metrics.csv")
    for name, counts in truth.items():
        g = P[P.plate == f"synthetic_{name}"]
        assert len(g) == fmt, f"{len(g)} of {fmt} wells analysed on {name}"
        assert abs(g.whole_well_equiv.sum() - sum(counts)) / sum(counts) < 0.2
    R = pd.read_csv(f"{out}/Results/statistics_contrasts.csv"); rr = R[R.metric == "Colony number"].estimate.iloc[0]
    assert 0.35 < rr < 0.7
