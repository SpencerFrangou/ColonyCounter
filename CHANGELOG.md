# Changelog

All notable changes to ColonyCounter are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/). A patch release fixes a bug without changing results, a minor
release adds a feature or changes a default in a way that can change results, a major release changes the
input format or the command line.

## [Unreleased]

### Added
- 12-, 24- and 96-well plates are detected and analysed; the plate geometry is selected from the
  `plate_format` setting of the plate map (`PLATE_FORMATS` in `cfa_pipeline.py`). Verified on synthetic
  plates; not yet checked on real scans of these formats.
- `make_synthetic_plates.py --format 6|12|24|96 --n-control N --n-treated N`.
- `pyproject.toml`: the repository is pip-installable and installs a `colonycounter` command
  (`template`, `check`, `run`, `calibrate`, `version`).
- `cfa_run.py --version`; every run writes `run_info.json` (software and dependency versions, settings,
  plate-map checksum) and prints the version in the report header.
- Unit tests for the plate map, the image-analysis helpers and the statistics; synthetic tests for the
  12- and 24-well formats; a GitHub Actions workflow running them on Linux, macOS and Windows with
  Python 3.10 to 3.13, and monthly against current dependency releases.
- `colonycounter/make_offline_wheelhouse.sh` for installing on computers without internet access.
- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue templates, `docs/statement_of_need.md`.

### Changed
- `.gitignore` restored (it was skipped by the browser upload of 1.0.0).
- `remove_small_objects` calls adapted to the scikit-image 0.26 API (`max_size`); behaviour unchanged,
  deprecation warning removed.
- The plate-map template writes the nominal growth area that matches the chosen plate format.

### Fixed
- The launchers copied `colonycounter/requirements.txt` into `code_used/`; the file is at the repository
  root.
- `USAGE.md` referred to a blank plate-map file and a project document that are not in the repository.

## [1.0.0] - 2026-10-29

First public release: well detection by rim fitting, optical-density segmentation with edge handling,
fillable plate map, design-aware statistics (quasi-Poisson rate ratios, log-scale fold changes, Holm and
Benjamini-Hochberg adjustment), publication figures, HTML report, one-click launchers for macOS and
Windows, synthetic example and end-to-end test. Archived at Zenodo, DOI 10.5281/zenodo.23780011.

[Unreleased]: https://github.com/SpencerFrangou/ColonyCounter/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/SpencerFrangou/ColonyCounter/releases/tag/v1.0.0
