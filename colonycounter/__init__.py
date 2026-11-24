"""ColonyCounter: quantification of colony formation (clonogenic) assays from flatbed scans.

The modules are written as scripts that can also be imported:
  cfa_pipeline   well detection, optical density, colony segmentation and per-well metrics
  cfa_platemap   fillable Excel plate map (template, reader, validation)
  cfa_run        driver: analysis of every plate, design-aware statistics, figures, HTML report
  cfa_calibrate  re-derivation of the OD threshold and rim fraction for a new scanner or stain
  cfa_analysis   figure style and statistics helpers
  cli            the `colonycounter` command (template, check, run, calibrate, version)
"""
__version__ = "1.1.0.dev0"
