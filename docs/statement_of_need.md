# Statement of need

The colony formation (clonogenic) assay measures the capacity of single cells to proliferate into colonies and is
the standard functional read-out for cell survival after a genetic perturbation, a drug or radiation. The read-out
is still, in most laboratories, a manual count of crystal violet stained colonies, or a count obtained with an
ImageJ macro that is written for one plate type and one scanner and is rarely shared. Manual counting is slow,
differs between people and days, ignores colony size and stain mass, and cannot use the edge of the well, where
meniscus and wall artefacts accumulate. Existing automated tools address parts of the problem: ColonyArea
(Guzmán et al. 2014) measures stained area in ImageJ but not colony number; OpenCFU and CellProfiler count
colonies in single images but leave well detection, edge handling, plate layout and statistics to the user;
CFAssay (Braselmann et al. 2015) provides the statistics but takes counts as input. None of them goes from a
folder of plate scans and a description of the experiment to tested, reportable numbers in one step.

ColonyCounter closes that gap for the most common set-up, a flatbed scan of a stained multi-well plate. It locates
each well by fitting the bright rim ring, converts the stain to optical density against a per-well flat field so
that integrated OD is proportional to dye mass, segments colonies with hysteresis thresholding and a watershed
split, and applies an explicit edge rule derived from the radial density of segmented objects (the outer 12% of the
radius is excluded; colonies are counted by centroid and sized only when fully inside). Every empirical constant
is documented with the measurement that produced it, and a calibration script re-derives the two that depend on
the scanner and stain.

The second half of the problem is that the statistics must match the layout. A fillable Excel plate map records
which cell line and condition each well carries and which condition is the control. From that map ColonyCounter
decides, for every comparison, whether plate can be entered as a blocking factor (both conditions on every plate),
only partly, or not at all (one condition per plate, in which case the confounding is stated in the report). Colony
number is modelled by a quasi-Poisson regression with the analysed area as offset, yielding a rate ratio with a
confidence interval; area and intensity by linear models on the log scale; primary endpoints are Holm-adjusted and
secondary endpoints carry Benjamini-Hochberg q-values. The report records every decision, so a reader can see not
only the p-value but the model behind it.

The intended users are cell biologists without programming experience: a double-click launcher creates the
environment, writes the plate map, runs the analysis and opens a report with QC overlays of every well, and the
same code runs from the command line or as a Python package for pipelines. The software was developed for a
breast cancer study at the University of Toledo and validated on 300 dpi scans of six 6-well plates from two cell
lines; synthetic plates with known colony numbers are included for testing, and the 12-, 24- and 96-well formats
are supported on that basis pending validation on real scans.

## References

Braselmann H, Michna A, Heß J, Unger K. CFAssay: statistical analysis of the colony formation assay. Radiat Oncol
2015;10:223. doi:10.1186/s13014-015-0529-y

Franken NAP, Rodermond HM, Stap J, Haveman J, van Bree C. Clonogenic assay of cells in vitro. Nat Protoc
2006;1:2315-2319. doi:10.1038/nprot.2006.339

Geissmann Q. OpenCFU, a new free and open-source software to count cell colonies and other circular objects. PLoS
One 2013;8:e54072. doi:10.1371/journal.pone.0054072

Guzmán C, Bagga M, Kaur A, Westermarck J, Abankwa D. ColonyArea: an ImageJ plugin to automatically quantify
colony formation in clonogenic assays. PLoS One 2014;9:e92444. doi:10.1371/journal.pone.0092444

Stirling DR, Swain-Bowden MJ, Lucas AM, Carpenter AE, Cimini BA, Goodman A. CellProfiler 4: improvements in speed,
utility and usability. BMC Bioinformatics 2021;22:433. doi:10.1186/s12859-021-04344-9
