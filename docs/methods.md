# ColonyCounter: methods and parameters

## Image analysis

**Well detection.** A crystal violet colour index, (R + B)/2 − G, thresholded by Otsu's method, gives a stained-field mask. Adjacent wells touch in this mask, so the well lattice (2 × 3, 3 × 4, 4 × 6 or 8 × 12) is placed from the mask extent and refined by centroid iteration; each well's bright rim ring is then located by ray casting (720 rays; brightest point between 0.85 and 1.12 R exceeding the interior median by 15 grey levels) and fitted by a robust circle fit with iterative trimming; all wells on a plate share the median radius. The pixel scale is read from the file DPI and checked against the well pitch (39.12 mm for 6-well plates); the two agreed within 2% on the validation scans.

**Optical density.** Crystal violet absorbs near 590 nm, so the green channel carries the signal. Each pixel is converted to OD = log10(I_background / I_green), where I_background is a per-well flat field from a 7 mm median filter. OD is additive in dye mass (Beer–Lambert), so integrated OD is the correct total-stain readout. Mean OD per colony is a stain-density measure and should not be compared across plates unless the design is blocked or a reference target is used.

**Segmentation.** Hysteresis thresholding on OD with the high threshold at 0.10 OD (seeds) and the low threshold at 0.06 OD (growth), hole filling, and a marker-controlled watershed on the distance transform (peak separation 0.4 mm) to split touching colonies. Objects below 0.15 mm equivalent diameter are discarded. At 300 dpi the practical detection floor is about 0.25 mm; at 600 dpi colonies at the 50-cell definition (about 0.14 mm) come within reach.

**Edge handling.** Only the disk r ≤ 0.88 R is analysed. A colony is counted if its centroid lies inside (an unbiased counting rule); its size and intensity are measured only if it lies entirely inside. Objects whose centroid lies in the 0.80 to 0.88 R band must have circularity ≥ 0.5. Whole-well colony number = counted colonies / analysed area × nominal growth area (950 mm² for a 6-well plate).

## How the two empirical parameters were derived, and how to re-derive them

*Edge exclusion.* With rim-fitted centres, the density of segmented objects per mm² is flat across the well interior and rises three- to eight-fold beyond about 0.88 R, where the meniscus arc, wall-attached cells seen obliquely and fragments of colonies climbing the wall are found; these objects have lower OD and lower circularity. On a CCD flatbed scanner the rims are also displaced laterally with distance from the lens axis (up to 2.8 mm on the validation scans), which is why the rim is fitted per well and why the plate should sit in the same place on the bed every time. ColonyArea (Guzmán et al. 2014) excludes 5% of the diameter for the same reason.

*OD threshold.* Per-well Otsu thresholds differ between cell lines with dense and faint colonies, so an automatic per-well threshold biases comparisons. A fixed threshold sweep (0.06 to 0.25 OD) shows a count plateau below about 0.11 OD on all validation plates and a steep fall above 0.14 OD for faintly staining lines; the OD noise floor of unstained pixels is below 0.02. A single absolute threshold of 0.10 OD was therefore fixed.

`python colonycounter/cfa_calibrate.py PlateMap.xlsx --out Calibration` writes the threshold sweep, the radial object-density profile and a plot of both for any plate map; choose the threshold at the low end of the count plateau and the rim fraction as the last radial band before the density rises, and enter both in the Settings sheet. Compare automated with manual counts of two or three wells once per cell line and stain protocol (Cai et al. 2011).

## Statistics

The well is the unit of replication. For each comparison the design is classified from the plate layout as blocked (plates carry both conditions; plate enters as a blocking factor), partly blocked, or confounded (one condition per plate; no plate term; confounding reported). Colony number: quasi-Poisson GLM with log(analysed area) offset and t-based inference on the estimated dispersion (rate ratios with 95% CI). Stained area fraction, integrated OD per mm², mean colony area: OLS on log values (fold changes). Mean diameter and mean colony OD: OLS (differences). Primary endpoints (colony number, stained area) are adjusted by Holm across all comparisons within a cell line; secondary endpoints by Benjamini–Hochberg (q-values). Games–Howell (three or more groups) or Welch t, Mann–Whitney, and Hedges' g with a bootstrap 95% CI are reported for every comparison; a Shapiro–Wilk test on model residuals flags comparisons where the rank-based p-value is the safer one. Groups with fewer than three wells receive estimates without p-values. Colony-level distributions are shown as descriptive summaries only, since colonies within a well are not independent.

Note that colony number measured at a single seeding density is a valid comparison at that density but cannot be converted into a clonogenic survival fraction when cells cooperate (Brix et al. 2020); a seeding-density series analysed with CFAcoop is needed for that.

## References

Franken NAP et al. Clonogenic assay of cells in vitro. Nat Protoc 2006. doi:10.1038/nprot.2006.339
Brix N et al. Analysis of clonogenic growth in vitro. Nat Protoc 2021. doi:10.1038/s41596-021-00615-0
Brix N et al. The clonogenic assay: robustness of plating efficiency-based analysis is strongly compromised by cellular cooperation. Radiat Oncol 2020. doi:10.1186/s13014-020-01697-y
Guzmán C et al. ColonyArea: an ImageJ plugin to automatically quantify colony formation in clonogenic assays. PLoS One 2014. doi:10.1371/journal.pone.0092444
Cai Z et al. Optimized digital counting colonies of clonogenic assays using ImageJ software and customized macros. Int J Radiat Biol 2011. doi:10.3109/09553002.2011.622033
Bewes JM et al. Automated cell colony counting and analysis using the circular Hough image transform algorithm. Phys Med Biol 2008. doi:10.1088/0031-9155/53/21/007
Brzozowska B et al. Freeware tool for analysing numbers and sizes of cell colonies. Radiat Environ Biophys 2019. doi:10.1007/s00411-018-00772-z
Braselmann H et al. CFAssay: statistical analysis of the colony formation assay. Radiat Oncol 2015. doi:10.1186/s13014-015-0529-y
Feoktistova M et al. Crystal violet assay for determining viability of cultured cells. Cold Spring Harb Protoc 2016. doi:10.1101/pdb.prot087379
