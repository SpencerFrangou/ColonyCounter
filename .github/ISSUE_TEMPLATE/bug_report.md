---
name: Bug report
about: Something gave a wrong result, failed, or behaved unexpectedly
labels: bug
---

**What happened**
A clear description of the problem. If a well or colony was mis-detected, attach the QC overlay
(`Results/QC_overlays/<plate>_overlay.jpg`), cropped if the scan cannot be shared.

**What you expected**

**How to reproduce**
The command or launcher used, the plate map settings (threshold OD, rim fraction, plate format), and if
possible a small image that triggers the problem (the synthetic generator in `examples/` can often be
adapted).

**Environment**
- ColonyCounter version (`colonycounter version`, or `run_info.json` in the results folder):
- Operating system and Python version:
- Scanner, resolution and file format of the scans:
