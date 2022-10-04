# targetedNOMe_2022

This repository contains all the code that was used to analyze Targeted NOMe data for our 2022 Nature Genetics Paper: [Long-range phasing of dynamic, tissue-specific and allele-specific regulatory elements]([https://www.nature.com/ng/](https://www.nature.com/articles/s41588-022-01188-8)).

The code is grouped by the general topic of analysis and will output the figures that were used in the paper. The scripts have dependencies such as necessary libraries or hard-coded paths to data files; these are declared at the beginning of each file.

The `pipeline/` folder contains the entire processing pipeline from FAST5 to RDATA.

The `util/` folder contains helper R functions that were commonly used in the analysis.
