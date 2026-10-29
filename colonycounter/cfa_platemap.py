"""Fillable plate map for colony formation assays.

    python cfa_platemap.py template <folder with scans> [--out PlateMap.xlsx] [--format 6]
    (fill in the workbook, then)
    python cfa_run.py PlateMap.xlsx --out Results

Workbook sheets
  PlateMap     one block per scanned plate: file name, cell line, and the well grid. Type the condition
               name into each well cell. Leave a well blank to ignore it. Write "Line:Condition" in a well
               to override the plate-level cell line (plates that carry two lines).
  Controls     one row per cell line: which condition is the control. If Comparisons is empty, every other
               condition of that line is compared with this control.
  Comparisons  optional explicit pairs: name, cell line, treatment, control. Takes precedence over Controls.
  Settings     analysis parameters (threshold OD, rim fraction, minimum colony diameter, plate format).
"""
import argparse, os, re, sys
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

FORMATS = {6: (2, 3), 12: (3, 4), 24: (4, 6), 96: (8, 12)}
IMG_EXT = (".jpg", ".jpeg", ".png", ".tif", ".tiff")
DEFAULT_SETTINGS = [("threshold_od", 0.10, "segmentation threshold (OD); 0.06 x 0.6 hysteresis low is derived from it"),
                    ("rim_fraction", 0.88, "fraction of the fitted rim radius analysed (outer band excluded)"),
                    ("min_diam_mm", 0.15, "smallest object kept (mm equivalent diameter)"),
                    ("plate_format", 6, "wells per plate: 6, 12, 24 or 96"),
                    ("growth_area_mm2", 950, "nominal growth area per well used for whole-well equivalents"),
                    ("all_pairs", "no", "yes = test every pair of conditions within a line, in addition to the comparisons listed")]

def make_template(folder, out, fmt, blank=0):
    rows, cols = FORMATS[fmt]
    if blank: files = [""] * blank
    else:
        files = sorted(f for f in os.listdir(folder) if f.lower().endswith(IMG_EXT))
        if not files: sys.exit(f"no images found in {folder}")
    wb = Workbook(); ws = wb.active; ws.title = "PlateMap"
    head = Font(bold=True); fill = PatternFill("solid", fgColor="E8EEF5"); well_fill = PatternFill("solid", fgColor="FFF8E1")
    thin = Side(style="thin", color="999999"); box = Border(left=thin, right=thin, top=thin, bottom=thin)
    ws["A1"] = "PLATE MAP. One block per scanned plate. File = exact image file name in this folder. Type the condition into each well; blank = ignore; 'Line:Condition' overrides the plate cell line. Add more blocks by copying one."; ws["A1"].font = head
    ws["A2"] = "Then fill Controls (control condition per cell line) and, if needed, Comparisons (explicit treatment vs control pairs). Settings holds the analysis parameters."
    r = 3
    for f in files:
        ws.cell(r, 1, "File").font = head; ws.cell(r, 2, f); ws.cell(r, 2).fill = (well_fill if blank else fill); ws.cell(r, 2).border = box
        ws.cell(r + 1, 1, "Cell line").font = head; ws.cell(r + 1, 2).fill = well_fill; ws.cell(r + 1, 2).border = box
        ws.cell(r + 2, 1, "Scan date").font = head; ws.cell(r + 2, 2).fill = well_fill; ws.cell(r + 2, 2).border = box
        ws.cell(r + 3, 1, "Cells seeded / well").font = head; ws.cell(r + 3, 2).fill = well_fill; ws.cell(r + 3, 2).border = box
        g = r + 5
        for c in range(cols): ws.cell(g, 2 + c, c + 1).font = head; ws.cell(g, 2 + c).alignment = Alignment(horizontal="center")
        for i in range(rows):
            ws.cell(g + 1 + i, 1, "ABCDEFGH"[i]).font = head
            for c in range(cols):
                cell = ws.cell(g + 1 + i, 2 + c); cell.fill = well_fill; cell.border = box; cell.alignment = Alignment(horizontal="center")
        r = g + rows + 3
    ws.column_dimensions["A"].width = 20
    for c in range(cols): ws.column_dimensions[get_column_letter(2 + c)].width = 22
    wc = wb.create_sheet("Controls"); wc.append(["Cell line", "Control condition"]); [setattr(x, "font", head) for x in wc[1]]
    wc.column_dimensions["A"].width = 20; wc.column_dimensions["B"].width = 24
    wp = wb.create_sheet("Comparisons"); wp.append(["Name", "Cell line", "Treatment", "Control"]); [setattr(x, "font", head) for x in wp[1]]
    for col in "ABCD": wp.column_dimensions[col].width = 24
    wsx = wb.create_sheet("Settings"); wsx.append(["Parameter", "Value", "Meaning"]); [setattr(x, "font", head) for x in wsx[1]]
    for row in DEFAULT_SETTINGS: wsx.append(list(row))
    wsx["B5"] = fmt; wsx.column_dimensions["A"].width = 18; wsx.column_dimensions["C"].width = 80
    wb.save(out); print(f"template written: {out}  ({len(files)} plates, {fmt}-well)")

def read_platemap(path):
    """Returns (wells DataFrame, controls dict, comparisons DataFrame, settings dict)."""
    wb = load_workbook(path, data_only=True)
    ws = wb["PlateMap"]; vals = [[c.value for c in row] for row in ws.iter_rows()]
    settings = {r[0]: r[1] for r in wb["Settings"].iter_rows(min_row=2, values_only=True) if r[0]}
    fmt = int(settings.get("plate_format", 6)); rows, cols = FORMATS[fmt]
    wells = []; i = 0; folder = os.path.dirname(os.path.abspath(path))
    while i < len(vals):
        if vals[i] and vals[i][0] == "File":
            fname = vals[i][1]
            if not fname or not str(fname).strip(): i += 1; continue
            fname = str(fname).strip(); line = vals[i + 1][1]; date = vals[i + 2][1]; seeded = vals[i + 3][1]; g = i + 5
            for ri in range(rows):
                for ci in range(cols):
                    v = vals[g + 1 + ri][1 + ci] if g + 1 + ri < len(vals) and 1 + ci < len(vals[g + 1 + ri]) else None
                    if v is None or str(v).strip() == "": continue
                    v = str(v).strip(); wl = line; cond = v
                    if ":" in v: wl, cond = [x.strip() for x in v.split(":", 1)]
                    wells.append(dict(file=fname, path=os.path.join(folder, fname), plate=os.path.splitext(fname)[0], well=f"{'ABCDEFGH'[ri]}{ci + 1}",
                                      row=ri, col=ci, line=(str(wl).strip() if wl else None), condition=cond, scan_date=date, cells_seeded=seeded))
            i = g + rows + 1
        else: i += 1
    W = pd.DataFrame(wells)
    controls = {str(r[0]).strip(): str(r[1]).strip() for r in wb["Controls"].iter_rows(min_row=2, values_only=True) if r[0] and r[1]}
    comps = pd.DataFrame([dict(name=r[0], line=str(r[1]).strip(), treatment=str(r[2]).strip(), control=str(r[3]).strip())
                          for r in wb["Comparisons"].iter_rows(min_row=2, values_only=True) if r[1] and r[2] and r[3]])
    return W, controls, comps, settings

def validate(W, controls, comps):
    problems = []
    if W.empty: problems.append("no wells labelled in PlateMap")
    if W.line.isna().any(): problems.append("some labelled wells have no cell line (fill 'Cell line' for the plate or use Line:Condition)")
    for p, g in W.groupby("plate"):
        if not os.path.exists(g.path.iloc[0]): problems.append(f"image not found: {g.path.iloc[0]}")
    dup = W.duplicated(["plate", "well"]);
    if dup.any(): problems.append("duplicate well entries: " + ", ".join(W[dup].plate + ":" + W[dup].well))
    for line, g in W.groupby("line"):
        conds = set(g.condition)
        if comps.empty or line not in set(comps.line):
            if line not in controls: problems.append(f"line {line}: no control condition in Controls and no Comparisons listed")
            elif controls[line] not in conds: problems.append(f"line {line}: control '{controls[line]}' is not a condition on any plate of this line")
        for c, gg in g.groupby("condition"):
            if len(gg) < 2: problems.append(f"line {line}, condition {c}: only {len(gg)} well (need at least 2; 3 or more for tests)")
    for r in comps.itertuples():
        sub = W[W.line == r.line]
        for c in (r.treatment, r.control):
            if c not in set(sub.condition): problems.append(f"comparison {r.name}: '{c}' is not a condition of line {r.line}")
    return problems

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    t = sp.add_parser("template"); t.add_argument("folder"); t.add_argument("--out", default=None); t.add_argument("--format", type=int, default=6)
    t.add_argument("--blank", type=int, default=0, help="write N empty plate blocks instead of one per image file")
    c = sp.add_parser("check"); c.add_argument("platemap")
    a = ap.parse_args()
    if a.cmd == "template": make_template(a.folder, a.out or os.path.join(a.folder, "PlateMap.xlsx"), a.format, a.blank)
    else:
        W, ctl, cp_, st = read_platemap(a.platemap); pr = validate(W, ctl, cp_)
        print(W.groupby(["line", "condition"]).size().rename("wells")); print("controls:", ctl); print("comparisons:", len(cp_))
        print("OK" if not pr else "PROBLEMS:\n  " + "\n  ".join(pr))
