# GitHub setup: ColonyCounter

Browser steps; sign in to GitHub as **Spencer** first. The repository contents are in this folder on
Frangou_2TB: `GitHub_Repos/ColonyCounter/` (a git repository with one commit) and, for drag-and-drop upload,
`GitHub_Repos/_UPLOAD_ColonyCounter_full/`. No research data is included anywhere: the only images are the
synthetic plates generated on demand by `examples/make_synthetic_plates.py`.

## 1. Create the repository
1. https://github.com/new -> Repository name: `ColonyCounter` -> Public.
2. Leave README / .gitignore / license unticked -> Create repository.

## 2. Upload (choose one)
**a. Drag and drop (no terminal).** On the new repo's quick-setup page click "uploading an existing file";
Finder -> `Frangou_2TB/GitHub_Repos/_UPLOAD_ColonyCounter_full`; select everything inside (Cmd+A) and drag it into
the browser; Commit changes. Folders (colonycounter/, docs/, examples/, tests/) are kept.

**b. Terminal, if git is installed on the Mac** (keeps the commit history):
```
cd "/Volumes/Frangou_2TB/GitHub_Repos/ColonyCounter"
git remote add origin https://github.com/SpencerFrangou/ColonyCounter.git
git branch -M main
git push -u origin main
```
(GitHub asks for the username and a personal access token the first time.)

## 3. About box (gear next to "About")
- Description: `Quantification of colony formation (clonogenic) assays from flatbed scans: well detection, OD-based colony segmentation, fillable plate map, design-aware statistics, publication figures.`
- Topics: `clonogenic-assay` `colony-formation` `crystal-violet` `image-analysis` `cancer-biology` `python` `opencv` `scikit-image` `reproducible-research`

## 4. Optional: citable DOI via Zenodo
https://zenodo.org -> Log in with GitHub -> GitHub tab -> switch ColonyCounter ON -> on GitHub make a release
(Releases -> Draft a new release -> tag `v1.0.0`). Zenodo mints a DOI; paste it into CITATION.cff and the README badge.

## 5. Verify
Open https://github.com/SpencerFrangou/ColonyCounter: the README should show the badges and the workflow figure,
and `python tests/test_synthetic.py` run from a fresh clone should print PASS.
