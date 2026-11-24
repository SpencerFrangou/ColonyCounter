"""Plate map: template -> fill -> read -> validate round trip, Line:Condition override, error detection."""
import os, sys, tempfile
import pytest
from openpyxl import load_workbook
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "colonycounter"))
from cfa_platemap import make_template, read_platemap, validate, FORMATS

def _fill(path, blocks):
    """blocks: {file name: (line, {well: text})}"""
    wb = load_workbook(path); ws = wb["PlateMap"]
    for row in ws.iter_rows():
        if row[0].value == "File" and row[1].value in blocks:
            line, wells = blocks[row[1].value]; r = row[0].row; ws.cell(r + 1, 2, line)
            for well, text in wells.items():
                ri = "ABCDEFGH".index(well[0]); ci = int(well[1:]) - 1; ws.cell(r + 6 + ri, 2 + ci, text)
    wb.save(path); return wb

def _touch(folder, *names):
    for n in names: open(os.path.join(folder, n), "wb").write(b"\xff\xd8\xff")   # a JPEG marker is enough for the listing

@pytest.mark.parametrize("fmt", sorted(FORMATS))
def test_template_has_one_block_per_image_and_format_settings(fmt):
    d = tempfile.mkdtemp(); _touch(d, "p1.jpg", "p2.tif"); out = os.path.join(d, "PlateMap.xlsx")
    make_template(d, out, fmt)
    wb = load_workbook(out); files = [r[1].value for r in wb["PlateMap"].iter_rows() if r[0].value == "File"]
    assert files == ["p1.jpg", "p2.tif"]
    settings = {r[0]: r[1] for r in wb["Settings"].iter_rows(min_row=2, values_only=True)}
    assert settings["plate_format"] == fmt and settings["growth_area_mm2"] == {6: 950, 12: 380, 24: 190, 96: 32}[fmt]

def test_blank_template_writes_empty_blocks():
    d = tempfile.mkdtemp(); out = os.path.join(d, "PlateMap.xlsx"); make_template(d, out, 6, blank=3)
    wb = load_workbook(out); assert sum(1 for r in wb["PlateMap"].iter_rows() if r[0].value == "File") == 3
    W, controls, comps, settings = read_platemap(out); assert W.empty   # blank File cells are skipped

def test_round_trip_controls_and_line_override():
    d = tempfile.mkdtemp(); _touch(d, "a.jpg", "b.jpg"); out = os.path.join(d, "PlateMap.xlsx"); make_template(d, out, 6)
    wb = _fill(out, {"a.jpg": ("MCF10A", {"A1": "Vector", "A2": "Vector", "A3": "Vector", "B1": "KO", "B2": "KO", "B3": "KO"}),
                     "b.jpg": ("MCF10A", {"A1": "KO", "A2": "KO", "A3": "KO", "B1": "SUM149:Vector", "B2": "SUM149:Vector", "B3": "SUM149:KO"})})
    wb["Controls"].append(["MCF10A", "Vector"]); wb["Controls"].append(["SUM149", "Vector"]); wb.save(out)
    W, controls, comps, settings = read_platemap(out)
    assert len(W) == 12 and set(W.line) == {"MCF10A", "SUM149"}
    assert (W[(W.plate == "b") & (W.well == "B1")].line == "SUM149").all() and (W[(W.plate == "b") & (W.well == "B1")].condition == "Vector").all()
    assert controls == {"MCF10A": "Vector", "SUM149": "Vector"} and comps.empty
    assert W.path.iloc[0] == os.path.join(d, "a.jpg")
    problems = validate(W, controls, comps)
    assert problems == [f"line SUM149, condition KO: only 1 well (need at least 2; 3 or more for tests)"]

def test_validate_reports_missing_control_and_missing_image():
    d = tempfile.mkdtemp(); _touch(d, "a.jpg"); out = os.path.join(d, "PlateMap.xlsx"); make_template(d, out, 6)
    _fill(out, {"a.jpg": ("L1", {"A1": "X", "A2": "X", "B1": "Y", "B2": "Y"})}); os.remove(os.path.join(d, "a.jpg"))
    W, controls, comps, settings = read_platemap(out); problems = validate(W, controls, comps)
    assert any("image not found" in p for p in problems) and any("no control condition" in p for p in problems)

def test_comparisons_sheet_takes_precedence():
    d = tempfile.mkdtemp(); _touch(d, "a.jpg"); out = os.path.join(d, "PlateMap.xlsx"); make_template(d, out, 6)
    wb = _fill(out, {"a.jpg": ("L1", {"A1": "ctl", "A2": "ctl", "A3": "ctl", "B1": "t1", "B2": "t2", "B3": "t2"})})
    wb["Comparisons"].append(["t2 vs ctl", "L1", "t2", "ctl"]); wb.save(out)
    W, controls, comps, settings = read_platemap(out)
    assert len(comps) == 1 and comps.treatment.iloc[0] == "t2"
    assert not any("no control condition" in p for p in validate(W, controls, comps))
