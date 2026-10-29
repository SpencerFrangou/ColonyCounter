"""Colony formation assay quantification from flatbed scans of crystal-violet stained plates.

Pipeline (parameters established on 300 dpi validation scans of crystal violet stained 6-well plates):
  1. Stained-field mask from the crystal-violet colour index (R+B)/2 - G, Otsu threshold.
  2. Plate orientation and well lattice (n_rows x n_cols) from the mask extent; lattice fit.
  3. Per-well rim-ring refinement: ray-cast the bright rim ring, robust circle fit, shared radius.
  4. Optical density OD = log10(I_bg / I_green) with per-well flat-field (61 px median).
  5. Segmentation: hysteresis threshold (high = T, low = 0.6 T) on OD, T either per-well Otsu
     or a fixed experiment-wide value; distance-transform watershed to split touching colonies.
  6. Edge handling: analysed disk r <= RIM * R (0.88); colonies counted by centroid inside the disk
     (unbiased counting); size/intensity measured on colonies fully inside the disk; objects with
     centroid in the 0.80-0.88 R band must have circularity >= 0.5.
  7. Per-colony and per-well metrics, QC overlay.
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, cv2, pandas as pd
from scipy import ndimage as ndi, optimize
from skimage import filters, morphology, measure, segmentation, feature

WELL_DIAM_MM = 34.8; WELL_PITCH_MM = 39.12; GROWTH_AREA_MM2 = 950.0
RIM = 0.88; BAND = 0.80; CIRC_MIN = 0.5; MIN_DIAM_MM = 0.15

def purple_index(im):
    b, g, r = cv2.split(im.astype(np.float32)); return (r + b) / 2 - g

def stained_mask(im):
    ps = cv2.GaussianBlur(purple_index(im), (0, 0), 2)
    thr = filters.threshold_otsu(ps[ps > 1])
    m = ps > thr
    m = morphology.remove_small_objects(m, 5000); m = ndi.binary_fill_holes(m); m = morphology.opening(m, morphology.disk(5))
    lab = measure.label(m); props = measure.regionprops(lab)
    big = [p.label for p in props if p.area > 0.25 * max(q.area for q in props)]
    return m & np.isin(lab, big), ps, thr

def lattice_wells(mask, n_wells=6):
    ys, xs = np.nonzero(mask)
    w = np.percentile(xs, 99.5) - np.percentile(xs, 0.5); h = np.percentile(ys, 99.5) - np.percentile(ys, 0.5)
    if n_wells == 6: n_rows, n_cols = (2, 3) if w > h else (3, 2)
    else: raise NotImplementedError
    extent = w if n_cols >= n_rows else h
    pitch = extent / ((max(n_rows, n_cols) - 1) + WELL_DIAM_MM / WELL_PITCH_MM)
    rad = pitch * WELL_DIAM_MM / WELL_PITCH_MM / 2
    cx0, cy0 = xs.mean(), ys.mean()
    H, W = mask.shape; yy, xx = np.mgrid[:H, :W]; wells = []
    for ri in range(n_rows):
        for ci in range(n_cols):
            cx = cx0 + (ci - (n_cols - 1) / 2) * pitch; cy = cy0 + (ri - (n_rows - 1) / 2) * pitch
            for _ in range(2):
                m = ((xx - cx) ** 2 + (yy - cy) ** 2 <= (0.9 * rad) ** 2) & mask
                if m.sum() == 0: break
                cy, cx = ndi.center_of_mass(m)
            wells.append(dict(name=f"{'ABC'[ri]}{ci + 1}", row=ri, col=ci, cx=float(cx), cy=float(cy), r=float(rad)))
    return wells, n_rows, n_cols, pitch

def ring_points(gray, cx, cy, R, n_ang=720, lo=0.85, hi=1.12, min_prom=15):
    pts = []; rs = np.arange(lo * R, hi * R, 0.5); rin = np.arange(0.5 * R, 0.8 * R)
    for a in np.linspace(0, 2 * np.pi, n_ang, endpoint=False):
        ca, sa = np.cos(a), np.sin(a)
        prof = ndi.map_coordinates(gray, [cy + rs * sa, cx + rs * ca], order=1, mode="nearest")
        ref = np.median(ndi.map_coordinates(gray, [cy + rin * sa, cx + rin * ca], order=1))
        k = int(np.argmax(prof))
        if prof[k] - ref > min_prom and 0 < k < len(rs) - 1: pts.append((cx + rs[k] * ca, cy + rs[k] * sa))
    return np.array(pts)

def fit_circle(pts, p0, r_fixed=None):
    f = (lambda p: np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1]) - p[2]) if r_fixed is None else \
        (lambda p: np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1]) - r_fixed)
    p = np.array(p0, float)
    for _ in range(4):
        p = optimize.least_squares(f, p, loss="soft_l1", f_scale=2.0).x
        e = np.abs(f(p)); keep = e < max(1.5, 2.5 * np.median(e)); pts = pts[keep]
    return p, pts

def image_dpi(path):
    try:
        from PIL import Image
        d = Image.open(path).info.get("dpi"); return float(d[0]) if d and d[0] > 50 else None
    except Exception: return None

def detect_wells(im, dpi=None):
    mask, ps, thr = stained_mask(im)
    wells, n_rows, n_cols, pitch = lattice_wells(mask)
    gray = cv2.GaussianBlur(cv2.cvtColor(im, cv2.COLOR_BGR2GRAY).astype(np.float32), (0, 0), 1.0)
    fits = []
    for w in wells:
        pts = ring_points(gray, w["cx"], w["cy"], w["r"])
        (cx, cy, r), kept = fit_circle(pts, [w["cx"], w["cy"], w["r"]])
        fits.append((cx, cy, r, pts))
    r_shared = float(np.median([f[2] for f in fits]))
    for w, (cx, cy, r, pts) in zip(wells, fits):
        (cx2, cy2), kept = fit_circle(pts, [cx, cy], r_fixed=r_shared)
        w.update(cx_lat=w["cx"], cy_lat=w["cy"], r_lat=w["r"], cx=float(cx2), cy=float(cy2), r=r_shared, r_free=float(r), n_ring_pts=len(kept))
    # pixel scale from the lattice pitch of the ring-fitted centres (less parallax than the mask extent)
    cs = np.array([[w["cx"], w["cy"]] for w in wells])
    d = [np.hypot(*(cs[i] - cs[j])) for i in range(6) for j in range(i + 1, 6)]
    pitch_px = np.min(d) if len(d) else pitch
    pitch_px = np.median([x for x in d if x < 1.3 * pitch_px])
    px_per_mm_pitch = pitch_px / WELL_PITCH_MM
    px_per_mm = dpi / 25.4 if dpi else px_per_mm_pitch
    return wells, dict(px_per_mm_pitch=px_per_mm_pitch, dpi=dpi, n_rows=n_rows, n_cols=n_cols, pitch_px=pitch_px, px_per_mm=px_per_mm, r_ring_px=r_shared, r_ring_mm=r_shared / px_per_mm)

def stain_od(sub, mask, px_per_mm=8.5):
    k = int(round(7.0 * px_per_mm)) | 1                     # 7 mm flat-field kernel, odd
    g8 = sub[:, :, 1]; bg = cv2.medianBlur(g8, k).astype(np.float32); g = g8.astype(np.float32) + 1.0
    bg = np.maximum(bg, g); od = np.log10(bg / g); od[~mask] = 0; return od, bg

def well_otsu(im, w, px_per_mm=8.5):
    """Otsu threshold of OD in the interior (r <= 0.85 R) of one well; used to set a global T."""
    sub, rr, x0, y0 = crop(im, w); od, _ = stain_od(sub, rr <= 1.0, px_per_mm)
    return float(filters.threshold_otsu(od[rr <= 0.85]))

def crop(im, w):
    R = w["r"]; pad = int(R * 1.05) + 5
    y0, y1 = int(w["cy"]) - pad, int(w["cy"]) + pad; x0, x1 = int(w["cx"]) - pad, int(w["cx"]) + pad
    sub = im[y0:y1, x0:x1]; ys, xs = np.mgrid[y0:y1, x0:x1]
    return sub, np.hypot(xs - w["cx"], ys - w["cy"]) / R, x0, y0

def analyse_well(im, w, px_per_mm, T=None, overlay=None, plate=""):
    sub, rr, x0, y0 = crop(im, w); R = w["r"]
    full = rr <= 1.0; disk = rr <= RIM
    od, bg = stain_od(sub, full, px_per_mm)
    t = float(filters.threshold_otsu(od[rr <= 0.85])) if T is None else T
    min_area_px = max(4, int((MIN_DIAM_MM * px_per_mm) ** 2 * np.pi / 4))
    fg = filters.apply_hysteresis_threshold(od, 0.6 * t, t) & full
    fg = morphology.remove_small_objects(fg, min_area_px // 2); fg = ndi.binary_fill_holes(fg)
    dist = ndi.distance_transform_edt(fg)
    pk = feature.peak_local_max(dist, min_distance=max(3, int(round(0.4 * px_per_mm))), labels=measure.label(fg), exclude_border=False)
    mk = np.zeros_like(fg, dtype=np.int32); mk[pk[:, 0], pk[:, 1]] = np.arange(1, len(pk) + 1)
    lab = segmentation.watershed(-dist, mk, mask=fg); lab = morphology.remove_small_objects(lab, min_area_px)
    cols = []; keep_lab = np.zeros(lab.max() + 1, bool); counted = 0
    ccx, ccy = w["cx"] - x0, w["cy"] - y0
    for p in measure.regionprops(lab, intensity_image=od):
        r_c = np.hypot(p.centroid[1] - ccx, p.centroid[0] - ccy) / R
        if r_c > RIM: continue
        circ = 4 * np.pi * p.area / max(p.perimeter, 1) ** 2
        if r_c >= BAND and circ < CIRC_MIN: continue
        r_max = np.max(np.hypot(p.coords[:, 1] - ccx, p.coords[:, 0] - ccy)) / R
        counted += 1; keep_lab[p.label] = True
        a_mm2 = p.area / px_per_mm ** 2
        cols.append(dict(plate=plate, well=w["name"], area_mm2=a_mm2, eq_diam_mm=2 * np.sqrt(a_mm2 / np.pi), mean_od=p.intensity_mean,
                         integrated_od=p.intensity_mean * p.area, circularity=circ, r_norm=r_c, fully_inside=bool(r_max <= RIM)))
    lab = np.where(keep_lab[lab], lab, 0)
    A = disk.sum() / px_per_mm ** 2
    C = pd.DataFrame(cols)
    Cin = C[C.fully_inside] if len(C) else C
    wd = dict(plate=plate, well=w["name"], n_colonies=counted, analysed_area_mm2=A, density_per_mm2=counted / A,
              whole_well_equiv=counted / A * GROWTH_AREA_MM2, area_fraction=float((lab > 0)[disk].mean()),
              total_area_mm2=float((lab > 0)[disk].sum() / px_per_mm ** 2), integrated_od=float(od[(lab > 0) & disk].sum()),
              integrated_od_per_mm2=float(od[(lab > 0) & disk].sum()) / A, threshold_od=t,
              mean_diam_mm=float(Cin.eq_diam_mm.mean()) if len(Cin) else np.nan, median_diam_mm=float(Cin.eq_diam_mm.median()) if len(Cin) else np.nan,
              mean_colony_area_mm2=float(Cin.area_mm2.mean()) if len(Cin) else np.nan, mean_colony_od=float(Cin.mean_od.mean()) if len(Cin) else np.nan,
              background_green=float(bg[disk].mean()), px_per_mm=px_per_mm)
    if overlay is not None:
        c = (int(w["cx"]), int(w["cy"]))
        cv2.circle(overlay, c, int(R), (255, 0, 0), 1); cv2.circle(overlay, c, int(R * RIM), (0, 0, 255), 2)
        cnts, _ = cv2.findContours((lab > 0).astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(overlay, [cc + np.array([x0, y0]) for cc in cnts], -1, (0, 200, 0), 1)
        cv2.putText(overlay, f"{w['name']} n={counted}", (c[0] - 55, int(c[1] - R - 8)), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 0, 255), 2)
    return C, wd

def analyse_plate(path, plate, T=None, out_overlay=None):
    im = cv2.imread(path); wells, info = detect_wells(im, dpi=image_dpi(path))
    ov = im.copy(); Cs, Ws = [], []
    for w in wells:
        C, wd = analyse_well(im, w, info["px_per_mm"], T=T, overlay=ov, plate=plate); Cs.append(C); Ws.append(wd)
    if out_overlay:
        xs = [w["cx"] for w in wells]; ys = [w["cy"] for w in wells]; R = info["r_ring_px"]
        y0, y1 = int(min(ys) - 1.3 * R), int(max(ys) + 1.3 * R); x0, x1 = int(min(xs) - 1.3 * R), int(max(xs) + 1.3 * R)
        cv2.imwrite(out_overlay, ov[max(y0, 0):y1, max(x0, 0):x1], [cv2.IMWRITE_JPEG_QUALITY, 90])
    return pd.concat(Cs, ignore_index=True), pd.DataFrame(Ws), pd.DataFrame(wells).assign(plate=plate), info
