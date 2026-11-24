# Contributing to ColonyCounter

Thank you for considering a contribution. Bug reports, questions about an unexpected result, requests for
a plate format or scanner that is not handled, and pull requests are all welcome.

## Reporting a problem

Open an issue at <https://github.com/SpencerFrangou/ColonyCounter/issues> using the bug report template.
The most useful report contains the ColonyCounter version (`colonycounter version`, or `run_info.json` from
the results folder), the operating system, what you ran, what happened and what you expected. A QC overlay
image from `Results/QC_overlays/` shows at once whether wells or colonies were mis-detected. Please do not
attach scans that you are not free to share; a cropped or synthetic image that reproduces the problem is
enough.

## Asking for help

Questions about using the software, choosing settings for a new scanner or stain, or interpreting the
statistics can also be opened as issues (label `question`). Expect a reply within a week.

## Proposing a change

1. Open an issue first for anything larger than a typo, so the design can be discussed before work starts.
2. Fork the repository and create a branch from `main`.
3. Install in development mode and run the tests:
   ```
   pip install -r requirements-dev.txt
   pip install -e .
   pytest tests/test_platemap.py tests/test_pipeline_units.py      # fast
   pytest tests/                                                   # includes the synthetic end-to-end runs, several minutes
   ```
4. Keep the analysis deterministic: no randomness without a fixed seed, no network access, no research
   data in the repository. The synthetic generator in `examples/` is the place for test images.
5. Add or update a test for the change, add a line under `Unreleased` in `CHANGELOG.md`, and update
   `USAGE.md` or `docs/methods.md` if the behaviour or a parameter changed.
6. Open a pull request. The GitHub Actions workflow runs the test suite on Linux, macOS and Windows.

A change that alters results for existing data (a default threshold, the edge rule, the statistical model)
is a minor version bump and must be described in the changelog with the reason and the expected size of
the effect.

## Code style

The code is written as compact scripts that are also importable. Follow the existing style of the file you
are editing, keep functions small, prefer explicit parameters over hidden state, and document any empirical
constant with how it was derived (see `docs/methods.md` for the two existing ones).

## Licence

By contributing you agree that your contribution is licensed under the MIT licence of this repository.
