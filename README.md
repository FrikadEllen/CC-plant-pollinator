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
* AllData_PlotLevel_BothYears.csv - All the processed plot-level datasets combined into 1 file, where values have been calculated or aggregated from the raw data to provide 1 value per plot per year. This is the data used in the plot-level regression analyses.
* AllData_MultiPlot_BothYears.csv - All the seed and nectar datasets combined into 1 file. This is the data used in the seed and nectar regression analyses. Note that due to the different sampling protocols used for the different datasets, there are many NA values that do not represent missing data (though there are also many true NAs in some datasets as well). 



## Experiment Methodology
Full details can be found in the published paper for this research.

We ran an outdoor, open-air, simulated climate-warming experiment on an arable farm in the UK (Stockbridge Technology Centre in North Yorkshire (53o49’N–1o9’W)) for two years, to investigate the impact of climate change on arable wildflowers and flower-visiting insects. The experiment consisted of 24 outdoor 2 x 2 m plots in an agricultural field, separated by 2m buffers, in a randomised block design with 6 replicates of four treatments: 1.5°C increase in temperature above ambient (‘Heat’); 40% increase in precipitation (‘Water’); warming and precipitation treatments combined (‘Heat+Water’); and ambient conditions (‘Control’). The heated plots were warmed with non-convective infrared heaters suspended 1.5 m above them operating continuously from the date of assembly (16/04/14 and 15/04/15) until end of sampling (19/08/14 and 18/08/2015). The plots and buffers were sown with spring wheat (Triticum aestivum cultivar Tybalt) and the plots were additionally sown with an arable wildflower seed mixture.

Sampling took place between the start of flowering in early June and the end of August (i.e. harvest) in 2014 and 2015. Seven sampling rounds were conducted in each year, with the dates matched as closely as possible to ensure even sampling between years. 

Multiple datasets were collected. 1) Floral resources: all flowering plant species were identified, and all floral units counted, in each plot during each sampling round. Additionally, nectar volumes were sampled from 3 species in 2015 (one sampling event per species). 2) Flower-visiting insects: each plot was observed for a total of 20 minutes per sampling round, during which, insect specimens feeding from flowers were captured using a hand-net and later identified to species. 3) Flower-visitor networks: a species interaction network was constructed for each plot, for each year, and network descriptors describing structure and complexity were calculated. 4) Wildflower seed set: seed heads of the 5 most abundant wildflower species were collected for analysis (one or two sampling events per species). Seeds were counted, a dry weight measurement of all seeds was taken, and average seed weight calculated (mg). 
