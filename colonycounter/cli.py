"""Command-line entry point installed as `colonycounter`.

    colonycounter template <folder> [--out PlateMap.xlsx] [--format 6|12|24|96] [--blank N]
    colonycounter check <PlateMap.xlsx>
    colonycounter run <PlateMap.xlsx> [--out Results]
    colonycounter calibrate <PlateMap.xlsx> [--out Calibration] [--threshold T]
    colonycounter version

Each subcommand runs the corresponding script (cfa_platemap.py, cfa_run.py, cfa_calibrate.py) with the same
arguments, so `python colonycounter/cfa_run.py ...` and `colonycounter run ...` behave identically.
"""
import os, runpy, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = {"template": "cfa_platemap.py", "check": "cfa_platemap.py", "run": "cfa_run.py", "calibrate": "cfa_calibrate.py"}

def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    from colonycounter import __version__
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__.strip()); return 0
    if argv[0] in ("version", "--version", "-V"):
        print(f"colonycounter {__version__}"); return 0
    cmd = argv[0]
    if cmd not in SCRIPTS:
        print(f"unknown command '{cmd}'. Commands: {', '.join(SCRIPTS)}, version", file=sys.stderr); return 2
    script = os.path.join(HERE, SCRIPTS[cmd])
    sys.argv = [script] + (argv if cmd in ("template", "check") else argv[1:])
    if HERE not in sys.path: sys.path.insert(0, HERE)
    try: runpy.run_path(script, run_name="__main__")
    except SystemExit as e: return int(e.code or 0) if not isinstance(e.code, str) else (print(e.code, file=sys.stderr) or 1)
    return 0

if __name__ == "__main__": sys.exit(main())
