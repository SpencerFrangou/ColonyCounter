"""Check that every package the analysis needs is installed; install whatever is missing.
Run automatically by the launchers before each analysis. Exit code 0 = ready, 1 = could not install."""
import importlib, subprocess, sys, os
REQ = {"numpy": "numpy>=1.26", "pandas": "pandas>=2.2", "scipy": "scipy>=1.13", "skimage": "scikit-image>=0.24", "cv2": "opencv-python-headless>=4.10",
       "statsmodels": "statsmodels>=0.14", "matplotlib": "matplotlib>=3.9", "openpyxl": "openpyxl>=3.1", "PIL": "Pillow>=10"}
missing = []
for mod, spec in REQ.items():
    try: importlib.import_module(mod)
    except Exception: missing.append(spec)
if not missing:
    print("all required packages present"); sys.exit(0)
print("missing packages:", ", ".join(missing))
wheels = os.path.join(os.path.dirname(os.path.abspath(__file__)), "wheels")
cmd = [sys.executable, "-m", "pip", "install", "-q"] + (["--no-index", "--find-links", wheels] if os.path.isdir(wheels) and os.listdir(wheels) else []) + missing
print("installing" + (" from offline cache" if "--no-index" in cmd else " from the internet") + " ...")
r = subprocess.run(cmd)
if r.returncode != 0:
    print("installation failed. Connect to the internet (or fill Code/wheels with make_offline_wheelhouse.sh) and run again."); sys.exit(1)
for mod in REQ:
    try: importlib.import_module(mod)
    except Exception as e: print(f"{mod} still cannot be imported: {e}"); sys.exit(1)
print("packages installed; ready"); sys.exit(0)
