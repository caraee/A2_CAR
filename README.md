# A2_CAR
Scripts for analysing data and generating figures for the A2-CAR-T cell manuscript. 
Figures 1, 3, and 5B, as well as supplementary figures S2-S7, were generated using the scripts found here. 

Each figure script contains code for wrangling the data (included as text files), 
performing the statistical analyses, and generating figures. Figures were then assembled using the patchwork package and saved as
svg files. Small changes were made to the svg files using Inkscape (https://inkscape.org/).

Figures 1 and S2-S3 use the raw data in the folder Figure1_data. Figures 3 and S4-S7 use the raw data in Figure3_data. 
Figure 5B uses the raw data in the folder Figure5B_data.

In all files, if values in the .Start and .End columns are not the same, columns ending in .Start are the lower limit
for unknown values and columns ending in .End are the upper limit. 147 mM was chosen as the upper limit for blood glucose
because it was the largest value anyone in the Kieffer lab had personally seen in samples. For C-peptide, the upper limit
is the lower limit of detection for the assay after dilution (in pg/mL).
