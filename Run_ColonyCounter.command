#!/bin/bash
# ColonyCounter: ONE double-click does everything.
#   - creates the private Python environment inside this folder if it does not exist yet
#   - checks the required packages; installs any that are missing (internet, or colonycounter/wheels offline cache)
#   - asks for the folder of scans; writes + opens PlateMap.xlsx on the first run for that folder,
#     analyses the plates and opens the report on the next run
cd "$(dirname "$0")"
echo "== ColonyCounter =="
if ! command -v python3 >/dev/null 2>&1; then
  echo "python3 was not found. Install Python 3.10 or later from https://www.python.org/downloads/ and run this again."; read -p "Press Enter to close"; exit 1; fi
if [ ! -d .venv ]; then echo "first run: creating the private environment (.venv)"; python3 -m venv .venv || { echo "could not create .venv"; read -p "Press Enter to close"; exit 1; }; fi
source .venv/bin/activate
python colonycounter/ensure_env.py || { read -p "Press Enter to close"; exit 1; }
FOLDER=$(osascript -e 'POSIX path of (choose folder with prompt "Choose the folder that contains the scanned plates")' 2>/dev/null)
[ -z "$FOLDER" ] && { echo "no folder chosen"; exit 0; }
FOLDER="${FOLDER%/}"; echo "folder: $FOLDER"
if [ ! -f "$FOLDER/PlateMap.xlsx" ]; then
  python colonycounter/cfa_platemap.py template "$FOLDER" && open "$FOLDER/PlateMap.xlsx"
  echo; echo "PlateMap.xlsx was created in that folder and opened. Fill in the wells, Controls (and Comparisons if needed), save it, then double-click Run_ColonyCounter.command again."
  read -p "Press Enter to close"; exit 0
fi
echo "== checking plate map =="; python colonycounter/cfa_platemap.py check "$FOLDER/PlateMap.xlsx" || { read -p "Fix the plate map and run again. Press Enter to close"; exit 1; }
STAMP=$(date +%Y-%m-%d_%H%M); OUT="$FOLDER/Results_$STAMP"; mkdir -p "$OUT"
echo "== running analysis -> $OUT =="
python colonycounter/cfa_run.py "$FOLDER/PlateMap.xlsx" --out "$OUT" 2>&1 | grep -v findfont | tee "$OUT/run_log.txt"
mkdir -p "$OUT/code_used" && cp colonycounter/*.py colonycounter/requirements.txt "$OUT/code_used/" && cp "$FOLDER/PlateMap.xlsx" "$OUT/PlateMap_used.xlsx"
open "$OUT/report.html"
echo; echo "Done. Results, figures, QC overlays, the plate map used and the code used are in: $OUT"
read -p "Press Enter to close"
