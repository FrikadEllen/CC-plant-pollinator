# CC-plant-pollinator
All data files and scripts for my simulated climate change field experiment on plant-pollinator interactions.


## Scripts
* CC-pp_Species-accumm-analysis.R - Species accummulation curves and calculation of Chao estimates of species richness, using package vegan.
* CC-pp_Community-analysis.r - PERMANOVA analysis of the plant and insect communities, using package vegan.
* CC-pp_Plot-level-analysis.R - Regression analyses of all plot-level variables (1 value per plot per year).
* CC-pp_Seeds-nectar-analysis.R - Regression analyses of all seed and nectar data (multiple values per plot per year).
* CC-pp_Figures.R - Code for creating the figures used in the current version of the manuscript.


## Data Files
* Raw_FlowerMatrix_BothYears.csv - Floral abundance values collected for each species in each plot, summed across the sample rounds. The first 24 rows are from 2014, 25-48 are from 2015.
* Raw_InteractionData_BothYears.csv - Raw flower-visitor interaction data collected during plot observations, formatted as a list of interactions (edge list).
* TreatmentYearData.txt - a simple key containing Year and Treatment columns, with rownames serving as Plot. Used in PERMANOVA to tell vegan which treatment each plot was assigned to, as the species data were in matrix format.
* AllData_PlotLevel_BothYears.csv - All the intermediate plot-level datasets combined into 1 file, where values have been calculated or aggregated from the raw data to provide 1 value per plot per year. This is the data used in the plot-level regression analyses.
* AllData_MultiPlot_BothYears.csv - All the seed and nectar datasets combined into 1 file. This is the data used in the seed and nectar regression analyses.Note that due to the different sampling protocols used for the different datasets, there are many NA values that do not represent missing data (though there are also many true NAs in some datasets as well). 
