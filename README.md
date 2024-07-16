# LongTerm_ClonalTracking
This code and data was used to generate the figures for the Long-Term Tracking of Hematopoietic Clonal Dynamics and Mutations in Nonhuman Primate Undergoing Transplantation of Lentivirally-Transduced Hematopoeitic Stem and Progenitor Cells BJH paper by the Dunbar Lab at the National Institutes of Health.

After updating directories, you should be able to use the pipeline after downloading the repository.

All results can be by running the barcodeanalysis.R script (main script) after properly configuring the loaders. Loaders are monkey specific, this is where you specify the counts and key (samples of interest) files. In the loader, also ensure that the counts matrix has the right number of rows per column. Imagine a matrix with the number of columns as cell types (T, B, Mono, Gr) and number of rows for each unique timepoint. If a cell type does not have a timepoint that another or multiple cell types have, then put an NA in that position.

Update the directories in both the main script and loaders.

For the unique clone number graphs in Figure 2 and autocorrelation graphs in Figure 3, there are figure specific scripts and data that do not feed into the main script. You will have to run these separately.

Source folder contains scripts performing calculations/graphing are referred to by the main script and loaders.

Results folder is empty, and this is where your outputs from the main script will be stored. Results.txt is just to have github have the folder. Can delete once repository is downloaded.

Supplemental Figure 2 data and code is provided. It is separate from the main script as well. Simply update the directories in the script for results.

For the bias plots in Supplemental Figure 7, all the animals should run fine aside from ZJ38 and ZH19. This is because there are missing timepoints for the pairwise comparisons between cell types of interest. Simply use the loaders specifically made for these plots and update the main script to use this loader instead of the default loaders for these two animals. Aside from ZJ38 and ZH19, the other animals have matching timepoints and should not yield an error. If for some reason it does, this is likely the cause. In the key file and loader, make sure that across all 4 cell types, the timepoints are shared. If it is not, then it is best to exclude that timepoint as the bias calculation cannot then be calculated.
