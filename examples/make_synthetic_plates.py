"""Generate synthetic plate scans so the pipeline can be tried without any real data.

    python examples/make_synthetic_plates.py --out examples/synthetic [--seed 0] [--format 6|12|24|96]
                                             [--n-control 60] [--n-treated 30]

Writes two plates (control and treated; the treated plate has half the colonies and smaller colonies by
default), as 300 dpi JPEG with DPI metadata, plus a filled PlateMap.xlsx and ground_truth.txt with the
number of colonies drawn in every well. Colonies are soft purple discs with log-normal sizes on a lavender
well floor with a bright rim ring, mimicking a crystal violet stained plate on a flatbed scanner. For
formats other than 6-well the default colony numbers are scaled to the well area. Intended for testing the
software, not for any biological inference.
"""
import argparse, os, sys, numpy as np
from PIL import Image
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "colonycounter"))
from cfa_pipeline import PLATE_FORMATS
DPI = 300; PX_MM = DPI / 25.4

def plate(fmt, n_per_well, diam_mm, seed, size=(2000, 1500)):
    f = PLATE_FORMATS[fmt]; rows, cols, pitch, R = f["rows"], f["cols"], f["pitch_mm"], f["well_diam_mm"] / 2
    rng = np.random.default_rng(seed); W, H = size
    img = np.full((H, W, 3), (235, 233, 236), np.float32)            # scanner bed
    yy, xx = np.mgrid[:H, :W].astype(np.float32); cx0, cy0 = W / 2, H / 2
    centres = [(cx0 + (c - (cols - 1) / 2) * pitch * PX_MM, cy0 + (r - (rows - 1) / 2) * pitch * PX_MM) for r in range(rows) for c in range(cols)]
    Rpx = R * PX_MM; counts = []
    for cx, cy in centres:
        d = np.hypot(xx - cx, yy - cy)
        img[d <= Rpx] = (208, 196, 232)                                       # lavender well floor
        ring = (d > Rpx - 2) & (d <= Rpx + 3); img[ring] = (250, 248, 252)     # bright rim ring
        wall = (d > Rpx + 3) & (d <= Rpx + min(12, 0.15 * Rpx)); img[wall] = (215, 212, 220)
        n = rng.poisson(n_per_well); counts.append(n)
        for _ in range(n):
            rr = 0.93 * Rpx * np.sqrt(rng.uniform(0, 1)); th = rng.uniform(0, 2 * np.pi); x, y = cx + rr * np.cos(th), cy + rr * np.sin(th)
            dia = rng.lognormal(np.log(diam_mm), 0.35) * PX_MM; s = dia / 2
            m = np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * (s / 1.6) ** 2)); m = np.clip(m * 1.4, 0, 1)
            dark = np.array((95, 60, 165), np.float32); img = img * (1 - m[..., None] * 0.95) + dark * m[..., None] * 0.95
    img += rng.normal(0, 2.5, img.shape); img = np.clip(img, 0, 255).astype(np.uint8)
    return Image.fromarray(img), counts

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="examples/synthetic"); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--format", type=int, default=6, choices=sorted(PLATE_FORMATS)); ap.add_argument("--n-control", type=int, default=None); ap.add_argument("--n-treated", type=int, default=None)
    a = ap.parse_args(); f = PLATE_FORMATS[a.format]
    scale = f["growth_area_mm2"] / PLATE_FORMATS[6]["growth_area_mm2"]
    n_ctl = a.n_control if a.n_control is not None else max(4, round(60 * scale)); n_trt = a.n_treated if a.n_treated is not None else max(2, round(30 * scale))
    diam = (0.9, 0.7) if f["well_diam_mm"] >= 15 else (0.5, 0.4)
    os.makedirs(a.out, exist_ok=True); truth = []
    for name, n, dia, s in [("control", n_ctl, diam[0], a.seed), ("treated", n_trt, diam[1], a.seed + 1)]:
        im, counts = plate(a.format, n, dia, s); im.save(f"{a.out}/synthetic_{name}.jpg", quality=92, dpi=(DPI, DPI)); truth.append((name, counts))
        print(f"synthetic_{name}.jpg ({a.format}-well): colonies drawn per well {counts}")
    from cfa_platemap import make_template
    from openpyxl import load_workbook
    make_template(a.out, f"{a.out}/PlateMap.xlsx", a.format)
    wb = load_workbook(f"{a.out}/PlateMap.xlsx"); ws = wb["PlateMap"]
    for row in ws.iter_rows():
        if row[0].value == "File":
            cond = "control" if "control" in row[1].value else "treated"; r = row[0].row; ws.cell(r + 1, 2, "SynthLine")
            for i in range(f["rows"]):
                for c in range(f["cols"]): ws.cell(r + 6 + i, 2 + c, cond)
    wb["Controls"].append(["SynthLine", "control"]); wb.save(f"{a.out}/PlateMap.xlsx")
    with open(f"{a.out}/ground_truth.txt", "w") as fh:
        for name, counts in truth: fh.write(f"{name}\t" + "\t".join(map(str, counts)) + "\n")
    print("plate map written:", f"{a.out}/PlateMap.xlsx")

if __name__ == "__main__": main()
