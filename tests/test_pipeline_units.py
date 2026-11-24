"""Unit checks of the image-analysis helpers and the design classification on small synthetic inputs."""
import os, sys
import numpy as np, pandas as pd, pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "colonycounter"))
import cfa_pipeline as cp
import cfa_run

def test_plate_formats_are_consistent():
    for n, f in cp.PLATE_FORMATS.items():
        assert f["rows"] * f["cols"] == n and f["well_diam_mm"] < f["pitch_mm"]
        assert abs(np.pi * (f["well_diam_mm"] / 2) ** 2 - f["growth_area_mm2"]) / f["growth_area_mm2"] < 0.12
    cp.set_plate_format(24); assert cp.N_WELLS == 24 and cp.GROWTH_AREA_MM2 == 190.0
    cp.set_plate_format(6); assert cp.WELL_PITCH_MM == 39.12

def test_fit_circle_recovers_centre_and_radius_with_outliers():
    rng = np.random.default_rng(0); th = rng.uniform(0, 2 * np.pi, 400)
    pts = np.c_[100 + 40 * np.cos(th), 80 + 40 * np.sin(th)] + rng.normal(0, 0.3, (400, 2))
    pts = np.vstack([pts, rng.uniform(0, 200, (40, 2))])                    # 10% outliers
    (cx, cy, r), kept = cp.fit_circle(pts, [95, 85, 35])
    assert abs(cx - 100) < 0.3 and abs(cy - 80) < 0.3 and abs(r - 40) < 0.3 and len(kept) < len(pts)
    (cx2, cy2), kept2 = cp.fit_circle(pts, [95, 85], r_fixed=40.0)
    assert abs(cx2 - 100) < 0.3 and abs(cy2 - 80) < 0.3

def test_stain_od_is_zero_on_flat_background_and_positive_on_stain():
    sub = np.full((201, 201, 3), 200, np.uint8); sub[90:110, 90:110, 1] = 100   # a dark square in the green channel
    mask = np.ones((201, 201), bool); od, bg = cp.stain_od(sub, mask, px_per_mm=8.0)
    assert np.allclose(od[:50, :50], 0, atol=1e-3)
    assert abs(od[100, 100] - np.log10(201 / 101)) < 0.02 and bg[100, 100] >= 199

def test_ring_points_find_a_bright_ring():
    yy, xx = np.mgrid[:300, :300].astype(np.float32); d = np.hypot(xx - 150, yy - 140)
    gray = np.full(d.shape, 120, np.float32); gray[(d > 59) & (d < 62)] = 250
    pts = cp.ring_points(gray, 148, 142, 58)
    assert len(pts) > 600 and abs(np.median(np.hypot(pts[:, 0] - 150, pts[:, 1] - 140)) - 60.5) < 1.5

def test_detect_wells_on_a_drawn_six_well_plate():
    cp.set_plate_format(6); px = 6.0; pitch = cp.WELL_PITCH_MM * px; R = cp.WELL_DIAM_MM / 2 * px
    H, W = 700, 900; img = np.full((H, W, 3), 235, np.uint8); yy, xx = np.mgrid[:H, :W]
    truth = [(W / 2 + (c - 1) * pitch, H / 2 + (r - 0.5) * pitch) for r in range(2) for c in range(3)]
    for cx, cy in truth:
        d = np.hypot(xx - cx, yy - cy); img[d <= R] = (232, 196, 208); img[(d > R - 2) & (d <= R + 2)] = 252
    wells, info = cp.detect_wells(img, dpi=px * 25.4)
    assert [w["name"] for w in wells] == ["A1", "A2", "A3", "B1", "B2", "B3"]
    for w, (cx, cy) in zip(wells, truth): assert abs(w["cx"] - cx) < 1.5 and abs(w["cy"] - cy) < 1.5
    assert abs(info["r_ring_px"] - R) < 3 and abs(info["px_per_mm_pitch"] - px) / px < 0.02

def _P(rows):
    return pd.DataFrame(rows, columns=["line", "condition", "plate"])

def test_classify_design():
    P = _P([("L", "c", "p1"), ("L", "t", "p1"), ("L", "c", "p2"), ("L", "t", "p2")]); assert cfa_run.classify_design(P, "L", "t", "c")[0] == "blocked"
    P = _P([("L", "c", "p1"), ("L", "t", "p2")]); assert cfa_run.classify_design(P, "L", "t", "c")[0] == "confounded"
    P = _P([("L", "c", "p1"), ("L", "t", "p1"), ("L", "t", "p2")]); assert cfa_run.classify_design(P, "L", "t", "c")[0] == "partly blocked"

def test_comparisons_for_controls_and_all_pairs():
    W = pd.DataFrame(dict(line=["L"] * 3, condition=["ctl", "a", "b"]))
    out = cfa_run.comparisons_for(W, {"L": "ctl"}, pd.DataFrame(), {})
    assert sorted(c["treatment"] for c in out) == ["a", "b"] and all(c["control"] == "ctl" for c in out)
    out = cfa_run.comparisons_for(W, {"L": "ctl"}, pd.DataFrame(), {"all_pairs": "yes"})
    assert len(out) == 3 and any({c["treatment"], c["control"]} == {"a", "b"} for c in out)
    comps = pd.DataFrame([dict(name="only b", line="L", treatment="b", control="ctl")])
    out = cfa_run.comparisons_for(W, {"L": "ctl"}, comps, {}); assert len(out) == 1 and out[0]["source"] == "Comparisons sheet"

def test_quasi_poisson_contrast_recovers_a_known_rate_ratio():
    rng = np.random.default_rng(1); n = 12
    d = pd.DataFrame(dict(cond=pd.Categorical(["c"] * n + ["t"] * n, categories=["c", "t"]), plate=["p"] * (2 * n), analysed_area_mm2=[700.0] * (2 * n)))
    d["n_colonies"] = np.r_[rng.poisson(80, n), rng.poisson(40, n)]
    est, lo, hi, p, model, m = cfa_run.fit_contrast(d, "n_colonies", "count", blocked=False)
    assert 0.4 < est < 0.62 and lo < 0.5 < hi and p < 1e-4 and model.startswith("quasi-Poisson")
