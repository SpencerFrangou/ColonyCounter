# ColonyCounter: how to use

ONE FILE DOES EVERYTHING
  macOS:   double-click  Run_ColonyCounter.command
  Windows: double-click  Run_ColonyCounter.bat
  Every time it runs it (1) creates the private Python environment inside this folder if it is not there yet,
  (2) checks that every required package is present and installs any that are missing (internet needed only
  for that; or the offline cache in colonycounter/wheels), then (3) asks for the folder of scans.
  Requirement: Python 3.10 or later installed (python.org). On a Mac the first double-click may need
  right-click > Open; if Terminal says "permission denied", run once:  chmod +x Run_ColonyCounter.command

FIRST RUN ON A FOLDER OF SCANS
  It writes PlateMap.xlsx into that folder and opens it. Fill it in and save:
    PlateMap sheet   one block per scan: cell line, date, cells seeded; type the condition into each well
                     (blank = ignore well; "Line:Condition" when a plate carries two cell lines)
    Controls sheet   the control condition for each cell line
    Comparisons      optional explicit treatment vs control pairs (replace the automatic ones for that line)
    Settings         threshold OD 0.10, rim fraction 0.88, min colony 0.15 mm, plate format, all_pairs
  A blank generic map with 12 empty plate blocks is in Examples/PlateMap_BLANK_template.xlsx.

SECOND RUN ON THE SAME FOLDER (and every run after that)
  Double-click Run_ColonyCounter.command again, choose the same folder. It checks the plate map, analyses every plate
  into Results_<date_time>/ inside that folder, and opens report.html. The Results folder also holds
  per_well_metrics.csv, colonies.csv, statistics_contrasts.csv, statistics_design.csv, Figure_<line>.png/svg/pdf,
  QC_overlays/, run_log.txt, the plate map used and the code used. Re-running never overwrites an earlier run.
  To re-analyse after editing the plate map, just run again.

PC VERSION
  Run_ColonyCounter.bat is the Windows version; the same project folder works on both systems. The .venv environment is
  specific to the computer that created it, so the first run on a new computer re-creates it (a few minutes).

RUNNING OFFLINE
  The analysis never uses the network. Only installing missing packages needs internet. For a computer that
  never has internet: on a connected computer with the same operating system and Python version run
  bash colonycounter/make_offline_wheelhouse.sh (fills colonycounter/wheels), copy the folder over; the launcher then installs
  from colonycounter/wheels automatically.

AFTER A CHANGE OF SCANNER, RESOLUTION, PLATE OR STAIN
  python colonycounter/cfa_calibrate.py <folder>/PlateMap.xlsx --out <folder>/Calibration
  then set threshold_od and rim_fraction in the Settings sheet from calibration.png (see PROJECT_DOCUMENTATION).

COMMAND LINE (no dialogs)
  source .venv/bin/activate            (Windows: .venv\Scripts\activate)
  python colonycounter/cfa_platemap.py template <folder> [--blank 12] [--format 6|12|24|96]
  python colonycounter/cfa_platemap.py check <folder>/PlateMap.xlsx
  python colonycounter/cfa_run.py <folder>/PlateMap.xlsx --out <folder>/Results
