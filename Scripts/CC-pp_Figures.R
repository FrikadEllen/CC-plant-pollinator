.libPaths(c("E:\\Rlib", .libPaths()))

library(ggplot2)
library(cowplot)
library(viridis)
library(reshape2)


#
######################### Functions ####

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  datac$Variable = paste0(measurevar)
  
  return(datac)
}

######




# 1. Data summarising
# 2. Flower data
# 3. Visitor data
# 4. Interaction data
# 5. Guild stacked bar
# 6. Combined plot-level plot
# 7. Seed data



######################### 1. Data summarising ####

# I will create some data files with all the variables summarised as means and standard errors, 
# including sample sizes. This could be useful for plotting, and will be useful for tables.
# make sure to load the summarySE function at the top of the script.


# Multiple values per plot data:
multi_df = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(multi_df)
# some of the columns have been classed as character, need to set all relevant ones to numeric:
multi_df[, c(3:4, 6:8, 10:22, 24:27)] = lapply(multi_df[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)

CFnectar_SE = summarySE(multi_df, measurevar="CF15_Nectar", groupvars=c("Treatment"), na.rm = TRUE)
CFnectar_SE
#    Treatment  N      mean        sd         se         ci    Variable
# 1    Control 12 0.2095016 0.1064228 0.03072161 0.06761781 CF15_Nectar
# 2       Heat 19 0.1806854 0.1279346 0.02935022 0.06166252 CF15_Nectar
# 3 Heat+Water 13 0.1616343 0.1202982 0.03336471 0.07269545 CF15_Nectar
# 4      Water 17 0.2834891 0.1537230 0.03728331 0.07903708 CF15_Nectar

multi_name_ls = list("CF14_Seeds", "CF14_AWeight", "CM14_Seeds", "CM14_AWeight", "RDN15_AWeight",
                     "SW15_Seeds", "SW15_AWeight", "CW15_Seeds", "CW15_AWeight", "CF15_Seeds", 
                     "CF15_AWeight", "CM15_Seeds", "CM15_AWeight", "CF15_Nectar", "RDN15_Nectar", "SW15_Nectar")

multi_SE = Reduce(rbind,
                  lapply(multi_name_ls, function(x) {
                    df = summarySE(multi_df, measurevar=x, groupvars=c("Treatment"), na.rm = TRUE)
                    df[, c(7, 1:6)]
                    }
                    )
                  )

head(multi_SE, n=10)
#        Variable  Treatment  N         mean           sd           se           ci
# 1    CF14_Seeds    Control 60 2.533333e+01 4.538598e+00 5.859305e-01 1.172444e+00
# 2    CF14_Seeds       Heat 60 2.121667e+01 4.274844e+00 5.518800e-01 1.104309e+00
# 3    CF14_Seeds Heat+Water 60 2.160000e+01 4.388738e+00 5.665836e-01 1.133731e+00
# 4    CF14_Seeds      Water 60 2.401667e+01 4.986394e+00 6.437407e-01 1.288122e+00
# 5  CF14_AWeight    Control 60 3.498969e-03 1.061605e-03 1.370526e-04 2.742416e-04
# 6  CF14_AWeight       Heat 60 3.741822e-03 8.868701e-04 1.144944e-04 2.291028e-04
# 7  CF14_AWeight Heat+Water 60 3.544514e-03 1.024918e-03 1.323163e-04 2.647643e-04
# 8  CF14_AWeight      Water 60 3.168955e-03 8.554803e-04 1.104420e-04 2.209940e-04
# 9    CM14_Seeds    Control 36 2.790833e+02 4.269585e+01 7.115975e+00 1.444620e+01
# 10   CM14_Seeds       Heat 36 2.306389e+02 7.376687e+01 1.229448e+01 2.495912e+01

# save this file:
write.csv(multi_SE, "Data/MultiPlot_Summary_Long.csv", row.names=FALSE)



# Version with the seed weight columns in mg instead of grams (so the values are easier to interpret)

multi_df$CF14_AWeight_mg = multi_df$CF14_AWeight *1000
multi_df$CM14_AWeight_mg = multi_df$CM14_AWeight *1000
multi_df$CF15_AWeight_mg = multi_df$CF15_AWeight *1000
multi_df$CM15_AWeight_mg = multi_df$CM15_AWeight *1000
multi_df$RDN15_AWeight_mg = multi_df$RDN15_AWeight *1000
multi_df$SW15_AWeight_mg = multi_df$SW15_AWeight *1000
multi_df$CW15_AWeight_mg = multi_df$CW15_AWeight *1000

multi_name_ls_mg = list("CF14_Seeds", "CF14_AWeight_mg", "CM14_Seeds", "CM14_AWeight_mg", "RDN15_AWeight_mg",
                     "SW15_Seeds", "SW15_AWeight_mg", "CW15_Seeds", "CW15_AWeight_mg", "CF15_Seeds", 
                     "CF15_AWeight_mg", "CM15_Seeds", "CM15_AWeight_mg", "CF15_Nectar", "RDN15_Nectar", "SW15_Nectar")

multi_SE_mg = Reduce(rbind,
                  lapply(multi_name_ls_mg, function(x) {
                    df = summarySE(multi_df, measurevar=x, groupvars=c("Treatment"), na.rm = TRUE)
                    df[, c(7, 1:6)]
                  }
                  )
)

head(multi_SE_mg, n=10)
#           Variable  Treatment  N       mean         sd         se         ci
# 1       CF14_Seeds    Control 60  25.333333  4.5385978  0.5859305  1.1724441
# 2       CF14_Seeds       Heat 60  21.216667  4.2748443  0.5518800  1.1043094
# 3       CF14_Seeds Heat+Water 60  21.600000  4.3887375  0.5665836  1.1337311
# 4       CF14_Seeds      Water 60  24.016667  4.9863939  0.6437407  1.2881221
# 5  CF14_AWeight_mg    Control 60   3.498969  1.0616048  0.1370526  0.2742416
# 6  CF14_AWeight_mg       Heat 60   3.741822  0.8868701  0.1144944  0.2291028
# 7  CF14_AWeight_mg Heat+Water 60   3.544514  1.0249177  0.1323163  0.2647643
# 8  CF14_AWeight_mg      Water 60   3.168955  0.8554803  0.1104420  0.2209940
# 9       CM14_Seeds    Control 36 279.083333 42.6958513  7.1159752 14.4461977
# 10      CM14_Seeds       Heat 36 230.638889 73.7668732 12.2944789 24.9591190

# save this file:
write.csv(multi_SE_mg, "Data/MultiPlot_Summary_mg_Long.csv", row.names=FALSE)



##

# Plot-level data:
plot_df = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(plot_df)
# some of the columns have been classed as character, need to set all relevant ones to numeric:
plot_df[, c(23, 29, 32)] = lapply(plot_df[, c(23, 29, 32)], as.numeric)

names(plot_df)
plot_name_ls = list("WConnectance", "GeneralityHL", "InteractionEvenness", "VulnerabilityLL", "Visits",
                    "Chao", "FRichness", "TotalFlowers", "VisitsOverFlowers")


plot_SE = Reduce(rbind,
                  lapply(plot_name_ls, function(x) {
                  df = summarySE(plot_df, measurevar=x, groupvars=c("Year", "Treatment"), na.rm = TRUE)
                  df[, c(8, 1:7)]
                  }
                  )
)

head(plot_SE, n=10)
#        Variable Year  Treatment N      mean         sd          se         ci
# 1  WConnectance    1    Control 6 0.1388407 0.02184846 0.008919596 0.02292855
# 2  WConnectance    1       Heat 6 0.1812097 0.04986176 0.020355978 0.05232671
# 3  WConnectance    1 Heat+Water 6 0.1820920 0.03452382 0.014094288 0.03623052
# 4  WConnectance    1      Water 6 0.1656831 0.01674574 0.006836419 0.01757358
# 5  WConnectance    2    Control 6 0.1962259 0.03758678 0.015344738 0.03944490
# 6  WConnectance    2       Heat 6 0.2292381 0.05842878 0.023853451 0.06131725
# 7  WConnectance    2 Heat+Water 6 0.2440621 0.02442599 0.009971867 0.02563350
# 8  WConnectance    2      Water 6 0.1843373 0.02832294 0.011562793 0.02972310
# 9  GeneralityHL    1    Control 6 1.5305220 0.23443033 0.095705780 0.24601954
# 10 GeneralityHL    1       Heat 6 1.5643391 0.28946246 0.118172556 0.30377223

# save this file:
write.csv(plot_SE, "Data/PlotLevel_Summary_Long.csv", row.names=FALSE)


# version with values for both years combined:
plot_SE_comb = Reduce(rbind,
                 lapply(plot_name_ls, function(x) {
                   df = summarySE(plot_df, measurevar=x, groupvars=c("Treatment"), na.rm = TRUE)
                   df[, c(7, 1:6)]
                 }
                 )
)
head(plot_SE_comb)
#       Variable  Treatment  N      mean         sd          se         ci
# 1 WConnectance    Control 12 0.1675333 0.04191961 0.012101149 0.02663445
# 2 WConnectance       Heat 12 0.2052239 0.05754110 0.016610686 0.03655987
# 3 WConnectance Heat+Water 12 0.2130770 0.04313140 0.012450964 0.02740439
# 4 WConnectance      Water 12 0.1750102 0.02422804 0.006994033 0.01539376
# 5 GeneralityHL    Control 12 1.4110251 0.25363967 0.073219466 0.16115496
# 6 GeneralityHL       Heat 12 1.4342823 0.28841603 0.083258537 0.18325080

# save this file:
write.csv(plot_SE_comb, "Data/PlotLevel_Summary_Combined-years_Long.csv", row.names=FALSE)



#
########################## 2. Flower data ####

# I need to plot floral abundance and nectar in a facetted plot, or use cowplot to put them together.
# I cannot find the old code for the thesis versions of these figures, so I need to recreate them from
# scratch.

# start with:
plot_SE = read.csv("Data/PlotLevel_Summary_Long.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
multi_SE = read.csv("Data/MultiPlot_Summary_Long.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
plot_df = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
multi_df = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")

multi_df[, c(3:4, 6:8, 10:22, 24:27)] = lapply(multi_df[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)
plot_df[, c(23, 29, 32)] = lapply(plot_df[, c(23, 29, 32)], as.numeric)

plot_df$Year = as.factor(plot_df$Year)
plot_SE$Year = as.factor(plot_SE$Year)


##### 2.1 Mean floral abundance "TotalFlowers" ####

ggplot(plot_SE[ which(plot_SE$Variable == "TotalFlowers"), ], 
       aes(x=Treatment, y=mean, fill=Year)) + 
  geom_bar(position=position_dodge(), stat="identity",
           colour="black" , size=.4) +      # Thinner lines=
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                size=.4,    # Thinner lines
                width=.3,
                position=position_dodge(.9)) +
  #geom_point() +
  xlab("Treatment") +
  ylab("Floral abundance") +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  #scale_y_continuous(expand = expand_scale(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0))
        #plot.margin=unit(c(2,5,2,2),"mm")




### Alternative using the full data and summary functions:
ggplot(plot_df, aes(x=Treatment, y=TotalFlowers, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=TotalFlowers), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0.4, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Floral abundance") +
  scale_fill_viridis_d(begin = 0.2, end =  0.8) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

  
ggplot(plot_df, aes(x=Treatment, y=TotalFlowers, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=TotalFlowers), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0.4, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Floral abundance") +
  scale_fill_manual(values = c("#a6dba0", "#af8dc3")) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))  

# other colour option (same darness): ("#a6dba0", "#c2a5cf")




#
##### 2.2 Nectar volumes ####

multi_df = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(multi_df)
multi_df[, c(3:4, 6:8, 10:22, 24:27)] = lapply(multi_df[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)

# I need to melt the nectar data into long format
nectar_long = melt(multi_df[, c(1, 2, 25:27)],
                    id.vars = c("Treatment", "Plot"),
                    variable.name = "Species", value.name = "Nectar")
head(multi_df[, c(1, 2, 25:27)])
head(nectar_long)
str(nectar_long)
unique(nectar_long$Species)
nectar_long$Species = as.character(nectar_long$Species)
nectar_long$Species = factor(nectar_long$Species, 
                                levels = c("CF15_Nectar", "RDN15_Nectar", "SW15_Nectar"),
                                labels = c("C. cyanus", "L. purpureum", "V. persica"))



ggplot(nectar_long, aes(x=Treatment, y=Nectar, fill=Species)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Nectar), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Floral nectar volume (\U003BCL)") +
  scale_fill_manual(values = c("#8da0cb", "#fc8d62", "#66c2a5"), labels = c("C. cyanus", "L. purpureum", "V. persica")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic",size = 9))) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0))



# See if facetting looks better:
ppi = 300
png("Figures/Nectar_facet.png", width=5*ppi, height=7*ppi, res=ppi)
#pdf("Figures/Nectar_facet.pdf", width=5, height=7)
ggplot(nectar_long, aes(x=Treatment, y=Nectar, fill=Species)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Nectar), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Floral nectar volume (\U003BCL)") +
  facet_wrap(~ Species, nrow =3, scales = "free_y") +
  scale_fill_manual(values = c("#8da0cb", "#fc8d62", "#66c2a5"), labels = c("C. cyanus", "L. purpureum", "V. persica")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic",size = 9))) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        strip.text = element_text(face = "italic"),
        legend.position = "bottom",
        legend.margin=margin(t=0, r=0, b=0, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))
dev.off()



####### Use viridis colours that match the seed plot
#### Extract the colours:
vir_lite = function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}
vir_lite(viridis(5), ds=0.6, dv=0.6)
# [1] "#A54CBB" "#899DD1" "#71D3CF" "#9EE9A2" "#FEF17C"

# re-lvel the factors to match the species ordering of the seed plots:
multi_df = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
multi_df[, c(3:4, 6:8, 10:22, 24:27)] = lapply(multi_df[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)
nectar_long = melt(multi_df[, c(1, 2, 25:27)],
                   id.vars = c("Treatment", "Plot"),
                   variable.name = "Species", value.name = "Nectar")
nectar_long$Species = as.character(nectar_long$Species)
nectar_long$Species = factor(nectar_long$Species, 
                             levels = c("CF15_Nectar", "SW15_Nectar", "RDN15_Nectar"),
                             labels = c("C. cyanus", "V. persica", "L. purpureum"))

# vertical facets:
ppi = 300
png("Figures/Nectar_facet_35-6.png", width=3.5*ppi, height=6*ppi, res=ppi)
#pdf("Figures/Nectar_facet_35-6.pdf", width=3.5, height=6)
ggplot(nectar_long, aes(x=Treatment, y=Nectar, fill=Species)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=0.4) +
  geom_point(aes(x = Treatment, y=Nectar), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.35) +
  xlab("Treatment") +
  ylab("Floral nectar volume (\U003BCL)") +
  facet_wrap(~ Species, nrow =3, scales = "free_y") +
  scale_fill_manual(values = c("#A54CBB", "#71D3CF", "#FEF17C"), labels = c("C. cyanus", "V. persica", "L. purpureum")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic",size = 9), nrow = 2)) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.title=element_text(size=10),
        legend.text = element_text(size=9.5),
        axis.text = element_text(colour = "black"),
        panel.spacing.y = unit(0.3, "lines"),
        axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 10, margin = margin(1,0,1,0, "mm")),
        axis.title.y = element_text(size=10, margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x = element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
        strip.text = element_text(face = "italic"),
        legend.position = "bottom",
        #legend.position = c(0.22, -0.13), 
        #legend.direction = "horizontal",
        legend.margin=margin(t=-2, r=0, b=0, l=-10),
        plot.margin=unit(c(1,1,1,1),"mm"))
dev.off()


# horizontal facets:
ppi = 300
#png("Figures/Nectar_facet_7-23.png", width=7*ppi, height=2.3*ppi, res=ppi)
png("Figures/Nectar_facet_7-24.png", width=7*ppi, height=2.4*ppi, res=ppi)
#png("Figures/Nectar_facet_7-25.png", width=7*ppi, height=2.5*ppi, res=ppi)
#pdf("Figures/Nectar_facet_7-25.pdf", width=7, height=2.5)
ggplot(nectar_long, aes(x=Treatment, y=Nectar, fill=Species)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=0.4) +
  geom_point(aes(x = Treatment, y=Nectar), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.4) +
  xlab("Treatment") +
  ylab("Floral nectar volume (\U003BCL)") +
  facet_wrap(~ Species, nrow =1, scales = "free_y") +
  scale_fill_manual(values = c("#A54CBB", "#71D3CF", "#FEF17C"), labels = c("C. cyanus", "V. persica", "L. purpureum")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic",size = 9))) +
  scale_x_discrete(labels = c("Control", "Heat", "Heat+\nWater", "Water")) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.title=element_text(size=10),
        legend.text = element_text(size=9.5),
        axis.text = element_text(colour = "black"),
       # panel.spacing.y = unit(0.3, "lines"),
        panel.spacing.x = unit(0.3, "lines"),
        axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 10, margin = margin(1,0,1,0, "mm")),
        axis.title.y = element_text(size=10, margin=margin(t=0, r=5, b=0, l=0)),
        axis.title.x = element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
        strip.text = element_text(face = "italic"),
        legend.position = "bottom",
        legend.margin=margin(t=-2, r=0, b=0, l=-10),
        plot.margin=unit(c(1,1,1,1),"mm"))
dev.off()



#
########################## 3. Visitor data ####


# I need to plot visitor abundance and visits per flower in a facetted plot, or use cowplot to put them together.
# I cannot find the old code for the thesis versions of these figures, so I need to recreate them from
# scratch.

# start with:
plot_df = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
plot_df[, c(23, 29, 32)] = lapply(plot_df[, c(23, 29, 32)], as.numeric)
plot_df$Year = as.factor(plot_df$Year)
str(plot_df)


##### 3.1 Visitor abundance ####

ggplot(plot_df, aes(x=Treatment, y=Visits, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Visits), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Flower-visitor abundance") +
  scale_fill_manual(values = c("#a6dba0", "#af8dc3")) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))  



#
##### 3.2 Visits per flower ####


ggplot(plot_df, aes(x=Treatment, y=VisitsOverFlowers, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=VisitsOverFlowers), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Frequency of visits to flowers") +
  scale_fill_manual(values = c("#a6dba0", "#af8dc3")) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm")) 


#
########################## 4. Interaction data ####


#  I want to see if it is worth plotting the significant interaction variables.

# start with:
plot_df = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
plot_df[, c(23, 29, 32)] = lapply(plot_df[, c(23, 29, 32)], as.numeric)
plot_df$Year = as.factor(plot_df$Year)
str(plot_df)


##### 4.1 Weighted connectance ####

ggplot(plot_df, aes(x=Treatment, y=WConnectance, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=WConnectance), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Weighted Connectance") +
  scale_fill_manual(values = c("#a6dba0", "#af8dc3")) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm")) 


#
##### 4.2 Interaction Evenness ####


ggplot(plot_df, aes(x=Treatment, y=InteractionEvenness, fill=Year)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=InteractionEvenness), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Interaction evenness") +
  scale_fill_manual(values = c("#a6dba0", "#af8dc3")) +
  #scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm")) 



#
########################## 5. Guild stacked bar ####

# I want to try making a stacked bar chart showing the different community composition between
# years, and possibly treatments (if it shows anything).


raw_int=read.csv("Data/Raw_InteractionData_BothYears.csv", header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
raw_int$Record_ID = 1:length(raw_int$Specimen_ID)
str(raw_int)
# 'data.frame':	3882 obs. of  8 variables:
#   $ Year       : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
# $ Plot       : int  2 2 16 1 11 23 23 9 9 15 ...
# $ Round      : int  1 1 1 1 1 1 1 1 1 1 ...
# $ Plant      : chr  "Veronica_persica" "Lamium_purpureum" "Stellaria_media" "Stellaria_media" ...
# $ Specimen_ID: chr  "Sphaerophoria_sp" "Bombus_pascuorum" "Eristalis_tenax" "Anthomyiidae_sp" ...
# $ Basic_Type : chr  "Hoverfly" "Bumblebee" "Hoverfly" "Fly" ...
# $ Treatment  : chr  "Heat+Water" "Heat+Water" "Control" "Control" ...
# $ Record_ID  : int  1 2 3 4 5 6 7 8 9 10 ...

aggregate(Record_ID ~ Basic_Type, data = raw_int, length)
# Basic_Type Record_ID
# 1      Bumblebee       792
# 2      Butterfly        68
# 3            Fly       413
# 4       Honeybee       932
# 5       Hoverfly      1632
# 6           Moth         1
# 7     Parasitoid         1
# 8    Social_Wasp         1
# 9   Solitary_Bee        41
# 10 Solitary_Wasp         1

# I am only interested in looking at broad groups of diptera and hymnenoptera.
target_types = c("Bumblebee", "Fly", "Honeybee", "Hoverfly", "Social_Wasp", "Solitary_Bee", "Solitary_Wasp")

raw_int_target = raw_int[which(raw_int$Basic_Type %in% target_types), ]
raw_int_target$Guild = ifelse(raw_int_target$Basic_Type == "Bumblebee", paste0("Wild_hym"), 
                              ifelse(raw_int_target$Basic_Type == "Solitary_Bee", paste0("Wild_hym"),
                                     ifelse(raw_int_target$Basic_Type == "Social_Wasp", paste0("Wild_hym"),
                                            ifelse(raw_int_target$Basic_Type == "Solitary_Wasp", paste0("Wild_hym"), 
                                                   raw_int_target$Basic_Type))))
aggregate(Record_ID ~ Basic_Type, data = raw_int_target, length)
#      Basic_Type Record_ID
# 1     Bumblebee       792
# 2           Fly       413
# 3      Honeybee       932
# 4      Hoverfly      1632
# 5   Social_Wasp         1
# 6  Solitary_Bee        41
# 7 Solitary_Wasp         1

aggregate(Record_ID ~ Guild, data = raw_int_target, length)
#      Guild Record_ID
# 1      Fly       413
# 2 Honeybee       932
# 3 Hoverfly      1632
# 4 Wild_hym       835


#
##### 5.1 By Year and Treatment #### 

guild_t_df = aggregate(Record_ID ~ Year + Treatment + Guild, data = raw_int_target, length)
# calculate proportion:
guild_t_df = transform(guild_t_df, Proportion = ave(Record_ID, Year, Treatment, FUN = function(x) x/sum(x)))
colnames(guild_t_df)[4] = "Abundance"
str(guild_t_df)
# 'data.frame':	32 obs. of  5 variables:
# $ Year      : int  2014 2015 2014 2015 2014 2015 2014 2015 2014 2015 ...
# $ Treatment : chr  "Control" "Control" "Heat" "Heat" ...
# $ Guild     : chr  "Fly" "Fly" "Fly" "Fly" ...
# $ Abundance : int  62 50 49 39 69 32 61 51 68 194 ...
# $ Proportion: num  0.1216 0.1 0.1384 0.0688 0.189 ...
head(guild_t_df)
#   Year  Treatment Guild Abundance Proportion
# 1 2014    Control   Fly        62 0.12156863
# 2 2015    Control   Fly        50 0.10000000
# 3 2014       Heat   Fly        49 0.13841808
# 4 2015       Heat   Fly        39 0.06878307
# 5 2014 Heat+Water   Fly        69 0.18904110
# 6 2015 Heat+Water   Fly        32 0.05860806

write.csv(guild_t_df, "Data/GuildNumbers_AllTreatments_BothYears.csv", row.names=FALSE)

guild_t_df$Year = as.factor(guild_t_df$Year)
guild_t_df$Guild = factor(guild_t_df$Guild, 
                          levels = c("Hoverfly","Fly","Honeybee","Wild_hym"),
                          labels = c("Syrphidae","Other Diptera","Apis mellifera","Other Hymenoptera"))


# Abdundace by year and treatment:
ggplot(guild_t_df, aes(x=Treatment, y=Abundance, fill=Guild)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlab("Treatment") +
  ylab("Visitor abundance") +
  scale_fill_manual(values=c("#008837", "#a6dba0", "#c2a5cf", "#7b3294")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Year, nrow = 2) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

# Proportion by year and treatment:
ggplot(guild_t_df, aes(x=Treatment, y=Proportion, fill=Guild)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlab("Treatment") +
  ylab("Visitor proportion") +
  scale_fill_manual(values=c("#F1605DFF", "#B63679FF", "#721F81FF", "#2D1160FF"), name = "Insect guild") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Year, nrow = 2) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))



#
##### 5.2 By Year #### 

guild_df = aggregate(Record_ID ~ Year + Guild, data = raw_int_target, length)
guild_df = transform(guild_df, Proportion = ave(Record_ID, Year, FUN = function(x) x/sum(x)))
colnames(guild_df)[3] = "Abundance"
str(guild_df)
# 'data.frame':	8 obs. of  4 variables:
# $ Year      : int  2014 2015 2014 2015 2014 2015 2014 2015
# $ Guild     : chr  "Fly" "Fly" "Honeybee" "Honeybee" ...
# $ Abundance : int  241 172 228 704 1015 617 191 644
# $ Proportion: num  0.1439 0.0805 0.1361 0.3294 0.606 ...
head(guild_df)
#   Year    Guild Abundance Proportion
# 1 2014      Fly       241 0.14388060
# 2 2015      Fly       172 0.08048666
# 3 2014 Honeybee       228 0.13611940
# 4 2015 Honeybee       704 0.32943379
# 5 2014 Hoverfly      1015 0.60597015
# 6 2015 Hoverfly       617 0.28872251

write.csv(guild_df, "Data/GuildNumbers_BothYears.csv", row.names=FALSE)

guild_df$Year = as.factor(guild_df$Year)
guild_df$Guild = factor(guild_df$Guild, 
                        levels = c("Hoverfly","Fly","Honeybee","Wild_hym"),
                        labels = c("Syrphidae","Other Diptera","Apis mellifera","Other Hymenoptera"))

# Abdundace by year:
ggplot(guild_df, aes(x=Year, y=Abundance, fill=Guild)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlab("Year") +
  ylab("Visitor abundance") +
  scale_fill_manual(values=c("#008837", "#a6dba0", "#c2a5cf", "#7b3294")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

# Proportion by year:
ggplot(guild_df, aes(x=Year, y=Proportion, fill=Guild)) + 
  geom_bar(position = "stack", stat = "identity") +
  xlab("Year") +
  ylab("Visitor proportion") +
  scale_fill_manual(values=c("#F1605DFF", "#B63679FF", "#721F81FF", "#2D1160FF"), name = "Insect guild") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=6, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=5, r=0, b=0, l=0)),
        legend.margin=margin(t=3, r=0, b=3, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))




#
########################## 6. Combined plot-level plot ####

plot_df = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(plot_df)
plot_df[, c(23, 29, 32)] = lapply(plot_df[, c(23, 29, 32)], as.numeric)
plot_df$Year = factor(plot_df$Year,
                      levels = c(1, 2),
                      labels = c("2014", "2015"))

guild_df = read.csv("Data/GuildNumbers_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(guild_df)
head(guild_df)
guild_df$Year = as.factor(guild_df$Year)
guild_df$Guild = factor(guild_df$Guild, 
                        levels = c("Hoverfly","Fly","Honeybee","Wild_hym"),
                        labels = c("Syrphidae","Other Diptera","Apis mellifera","Other Hymenoptera"))



# I want to see if I can use cowplot to combine some of the figures.

flowers_p = ggplot(plot_df, aes(x=Treatment, y=TotalFlowers, fill=Year)) +
  geom_bar(position=position_dodge(0.8), width = 0.8, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=TotalFlowers), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.8), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Floral abundance") +
  scale_fill_manual(values = c("#FEAF77FF", "#FCFDBFFF")) + # lighted 2 magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))

visits_p = ggplot(plot_df, aes(x=Treatment, y=Visits, fill=Year)) +
  geom_bar(position=position_dodge(0.8), width = 0.8, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Visits), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.8), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Insect abundance") +
  scale_fill_manual(values = c("#FEAF77FF", "#FCFDBFFF")) + # lighted 2 magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))

visits_flowers_p = ggplot(plot_df, aes(x=Treatment, y=VisitsOverFlowers, fill=Year)) +
  geom_bar(position=position_dodge(0.8), width = 0.8, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=VisitsOverFlowers), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.8), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Visits per flower") +
  scale_fill_manual(values = c("#FEAF77FF", "#FCFDBFFF")) + # lighted 2 magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))

wconn_p = ggplot(plot_df, aes(x=Treatment, y=WConnectance, fill=Year)) +
  geom_bar(position=position_dodge(0.8), width = 0.8, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=WConnectance), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.8), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Weighted\nconnectance") +
  scale_fill_manual(values = c("#FEAF77FF", "#FCFDBFFF")) + # lighted 2 magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))

inteven_p = ggplot(plot_df, aes(x=Treatment, y=InteractionEvenness, fill=Year)) +
  geom_bar(position=position_dodge(0.8), width = 0.8, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=InteractionEvenness), shape = 21, size = 1.1,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.8)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.8), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Interaction evenness") +
  scale_fill_manual(values = c("#FEAF77FF", "#FCFDBFFF")) + # lighted 2 magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))

guild_p = ggplot(guild_df, aes(x=Year, y=Proportion, fill=Guild)) + 
  geom_bar(position = "stack", stat = "identity", width = 0.8) +
  xlab("Year") +
  ylab("Insect proportion") +
  scale_fill_manual(values=c("#F1605DFF", "#B63679FF", "#721F81FF", "#2D1160FF"), name = "Insect guild") + # darkest 4 reverse magma(7)
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.0)), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.margin=margin(t=3, r=0, b=3, l=0))


plot_panel = plot_grid(flowers_p + theme(legend.position="none",
                                         axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                         axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                         axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
                                         plot.margin=unit(c(4,2,0,3),"mm"), # left side
                                         panel.border = element_blank(),
                                         axis.line    = element_line(color='black'),
                                         axis.text = element_text(colour = "black"),
                                         panel.grid.major.x = element_blank()), 
                       visits_p + theme(legend.position="none",
                                        axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                        axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                        axis.title.y=element_text(size=10, margin=margin(t=0, r=3, b=0, l=0)),
                                        plot.margin=unit(c(4,2,0,3),"mm"), # right
                                        panel.border = element_blank(),
                                        axis.line    = element_line(color='black'),
                                        axis.text = element_text(colour = "black"),
                                        panel.grid.major.x = element_blank()), 
                       visits_flowers_p + theme(legend.position="none",
                                                axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                                axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                                axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
                                                plot.margin=unit(c(4,2,0,3),"mm"), # left side
                                                panel.border = element_blank(),
                                                axis.line    = element_line(color='black'),
                                                axis.text = element_text(colour = "black"),
                                                panel.grid.major.x = element_blank()), 
                       guild_p + theme(legend.title=element_text(size=10),
                                       legend.text = element_text(size=9.5),
                                       axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                       axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                       axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
                                       plot.margin=unit(c(4,1,0,3),"mm"), # right
                                       panel.border = element_blank(),
                                       axis.line    = element_line(color='black'),
                                       axis.text = element_text(colour = "black"),
                                       panel.grid.major.x = element_blank()),
                       wconn_p + theme(legend.position="none",
                                       axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                       axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                       axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
                                       plot.margin=unit(c(4,2,0,3),"mm"), # left side
                                       panel.border = element_blank(),
                                       axis.line    = element_line(color='black'),
                                       axis.text = element_text(colour = "black"),
                                       panel.grid.major.x = element_blank()), 
                       inteven_p + theme(legend.position="none",
                                         axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
                                         axis.title.x=element_text(size=10, margin=margin(t=5, r=0, b=0, l=0)),
                                         axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
                                         plot.margin=unit(c(4,2,0,3),"mm"), # right
                                         panel.border = element_blank(),
                                         axis.line    = element_line(color='black'),
                                         axis.text = element_text(colour = "black"),
                                         panel.grid.major.x = element_blank()),
                       #labels = "AUTO", 
                       labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
                       #label_fontface = "italic",
                       label_size = 11, #scale = 1, label_colour = "blue",
                       #label_x = 0, label_y = 0.11,
                       #vjust = -0.005,
                       align = "v", axis = 'l',
                       ncol = 2)

panel_legend = get_legend(
  flowers_p + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.margin=margin(t=-1, r=0, b=0, l=0),
          legend.title=element_text(size=10),
          plot.margin=unit(c(0,0,0,0),"mm"))
)


ppi = 300
png("Figures/Plot-level_cowplot_7-6.png", width=7*ppi, height=6*ppi, res=ppi)
#pdf("Figures/Plot-level_cowplot_7-6.pdf", width=7, height=6)
plot_grid(plot_panel, panel_legend, ncol = 1, rel_heights = c(1, 0.05))
dev.off()



#
########################## 7. Seed data ####


# I want to see if it is worth plotting the significant seed variables in a big facet plot.

multi_df = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(multi_df)
# some of the columns have been classed as character, need to set all relevant ones to numeric:
multi_df[, c(3:4, 6:8, 10:22, 24:27)] = lapply(multi_df[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)

# I need to  melt the seed data into long format dfs, one for seed number, one for seed weight
str(multi_df[, c(1, 2, 3, 7, 15, 17, 19, 21)])
seed_n_long = melt(multi_df[, c(1, 2, 3, 7, 15, 17, 19, 21)],
                   id.vars = c("Treatment", "Plot"),
                   variable.name = "Species", value.name = "Seed_Number")

str(multi_df[, c(1, 2, 4, 8, 14, 16, 18, 20, 22)])
seed_w_long = melt(multi_df[, c(1, 2, 4, 8, 14, 16, 18, 20, 22)],
                   id.vars = c("Treatment", "Plot"),
                   variable.name = "Species", value.name = "Seed_Weight")
seed_w_long$Seed_W_mg = seed_w_long$Seed_Weight * 1000
head(seed_w_long)


#
#### 7.1 Seed number ####

# I need to relevel the species factor, and create a new factor for colouring (with fewer levels)
seed_n_long$Species = as.character(seed_n_long$Species)
unique(seed_n_long$Species)
# [1] "CF14_Seeds" "CM14_Seeds" "SW15_Seeds" "CW15_Seeds" "CF15_Seeds" "CM15_Seeds"
seed_n_long$Species_colour = seed_n_long$Species
seed_n_long$Species_colour = factor(seed_n_long$Species_colour, 
                                    levels = c("CF14_Seeds", "CF15_Seeds", "CM14_Seeds", "CM15_Seeds", "SW15_Seeds", "CW15_Seeds"),
                                    labels = c("C. cyanus", "C. cyanus", "G. segetum", "G. segetum", "V. persica", "S. media"))


seed_n_long$Species = factor(seed_n_long$Species, 
                             levels = c("CF14_Seeds", "CF15_Seeds", "CM14_Seeds", "CM15_Seeds", "SW15_Seeds", "CW15_Seeds"),
                             labels = c("C. cyanus 2014", "C. cyanus 2015", "G. segetum 2014", "G. segetum 2015", "V. persica", "S. media"))

# this makes the species names itlaic, but leaves the year not italic:
levels(seed_n_long$Species) = c("italic('C. cyanus ')*(2014)", "italic('C. cyanus ')*(2015)",
                                "italic('G. segetum ')*(2014)", "italic('G. segetum ')*(2015)",
                                "italic('V. persica')", "italic('S. media')")


# plot the seed number in a facet:
ggplot(seed_n_long, aes(x=Treatment, y=Seed_Number, fill=Species_colour)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Seed_Number), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Seeds per seedhead") +
  facet_wrap(~ Species, ncol =3, scales = "free_y", labeller=label_parsed) + # need the labeller function to parse the italic text
  #scale_fill_manual(values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF"), name = "Species") + # viridis colours
  scale_fill_manual(values = c("#0072B2", "#F0E442", "#009E73", "#999999"), name = "Species") + # okabe colours
  guides(fill = guide_legend(label.theme = element_text(face = "italic", size = 9))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=4, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=4, r=0, b=0, l=0)),
        #legend.position = "bottom",
        legend.margin=margin(t=0, r=0, b=0, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

############# I need to tidy this plot up:
# change the colours
# sort out the theme (axis text black and bigger, axis title smaller?)
# then make the equivalent plot for seed weight

# then use cowplot to put the two plots together. Put the legend in the empty space on the seed number plot, as it has 1
# fewer panels (weight has RDN). This might need some fiddling to get right (or it might be easy).

# viridis(5) colours:
# "#440154FF" "#3B528BFF" "#21908CFF" "#5DC863FF" "#FDE725FF"
# viridis(6) colours:
# "#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"

# okabe colours:
# cornflower "#0072B2"
# corn marigold "#F0E442" (pale yellow) or "#E69F00" (orange)
# s media "#999999" (grey) or "#009E73" (green) or "#F0E442" (pale yellow)
# v persica "#56B4E9" (pale blue) or "#009E73" (green)
# l purp "#D55E00" (dark orange) or "#CC79A7" (pink)


#
#### 7.2 Seed weight ####

# I need to relevel the species factor, and create a new factor for colouring (with fewer levels)
seed_w_long$Species = as.character(seed_w_long$Species)
unique(seed_w_long$Species)
# [1] "CF14_AWeight"  "CM14_AWeight"  "RDN15_AWeight" "SW15_AWeight"  "CW15_AWeight"  "CF15_AWeight"  "CM15_AWeight" 
seed_w_long$Species_colour = seed_w_long$Species
seed_w_long$Species_colour = factor(seed_w_long$Species_colour, 
                                    levels = c("CF14_AWeight", "CF15_AWeight", "CM14_AWeight", "CM15_AWeight", 
                                               "SW15_AWeight", "CW15_AWeight", "RDN15_AWeight"),
                                    labels = c("C. cyanus", "C. cyanus", "G. segetum", "G. segetum", 
                                               "V. persica", "S. media", "L. purpureum"))

seed_w_long$Species = factor(seed_w_long$Species, 
                             levels = c("CF14_AWeight", "CF15_AWeight", "CM14_AWeight", "CM15_AWeight", 
                                        "SW15_AWeight", "CW15_AWeight", "RDN15_AWeight"),
                             labels = c("C. cyanus 2014", "C. cyanus 2015", "G. segetum 2014", "G. segetum 2015", 
                                        "V. persica", "S. media", "L. purpureum"))

# this makes the species names itlaic, but leaves the year not italic:
levels(seed_w_long$Species) = c("italic('C. cyanus ')*(2014)", "italic('C. cyanus ')*(2015)",
                                "italic('G. segetum ')*(2014)", "italic('G. segetum ')*(2015)",
                                "italic('V. persica')", "italic('S. media')", "italic('L. purpureum')")




# plot seed weight in a facet:
sw_plot = ggplot(seed_w_long, aes(x=Treatment, y=Seed_Weight, fill=Species_colour)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean,
           colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Seed_Weight), shape = 21, 
             position = position_jitterdodge(jitter.width = 0.2, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), 
                size=0.4, width = 0.3) +
  xlab("Treatment") +
  ylab("Seed weight") +
  facet_wrap(~ Species, ncol =3, scales = "free_y", labeller=label_parsed) + # need the labeller function to parse the italic text
  #scale_fill_manual(values = c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF"), name = "Species") + # viridis colours
  scale_fill_manual(values = c("#0072B2", "#F0E442", "#009E73", "#999999", "#CC79A7"), name = "Species") + # okabe colours
  guides(fill = guide_legend(label.theme = element_text(face = "italic", size = 9), ncol=2)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_bw() +
  theme(axis.title.y = element_text(margin=margin(t=0, r=4, b=0, l=0)),
        axis.title.x = element_text(margin=margin(t=4, r=0, b=0, l=0)),
        #legend.position = "bottom",
        legend.margin=margin(t=0, r=0, b=0, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

lemon::reposition_legend(sw_plot, 'top left', panel='panel-3-3')

lemon::reposition_legend(sw_plot, 'center', panel=c('panel-3-2','panel-3-3'))


#
#### 7.3 Combine Number and weight using cowplot ####

vir_lite = function(cols, ds=0.4, dv=0.7) {
  cols = rgb2hsv(col2rgb(cols))
  cols["v", ] = cols["v", ] + dv*(1 - cols["v", ])
  cols["s", ] = ds*cols["s", ]
  apply(cols, 2, function(x) hsv(x[1], x[2], x[3]))
}

vir_lite(viridis(5), ds=0.7, dv=0.7)

sn_plot = ggplot(seed_n_long, aes(x=Treatment, y=Seed_Number, fill=Species_colour)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=0.4) +
  geom_point(aes(x = Treatment, y=Seed_Number), shape = 21, size = 1,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Seeds per seedhead") +
  facet_wrap(~ Species, ncol =3, scales = "free_y", labeller=label_parsed) + # need the labeller function to parse the italic text
  scale_fill_manual(values = vir_lite(viridis(5), ds=0.6, dv=0.6), name = "Species") +
  scale_x_discrete(labels = c("Control", "Heat", "Heat+\nWater", "Water")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic", size = 9))) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() + 
  theme(legend.position="none",
        axis.text = element_text(colour = "black"),
        panel.spacing.y = unit(0.2, "lines"),
        panel.spacing.x = unit(0.3, "lines"),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 10, margin = margin(0,0,0,0, "mm")),
        axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
        axis.title.x=element_text(size=10, margin=margin(t=4, r=0, b=0, l=0)),
        axis.title.y=element_text(size=10, margin=margin(t=0, r=4, b=0, l=0)),
        plot.margin=unit(c(1,1,1,1),"mm"))

sw_plot = ggplot(seed_w_long, aes(x=Treatment, y=Seed_W_mg, fill=Species_colour)) +
  geom_bar(position=position_dodge(0.7), width = 0.7, stat = 'summary', fun = mean, colour="black" , size=.4) +
  geom_point(aes(x = Treatment, y=Seed_W_mg), shape = 21, size = 1,
             position = position_jitterdodge(jitter.width = 0.3, jitter.height=0, dodge.width=0.7)) +
  geom_errorbar(stat = 'summary', fun.data = mean_se, position=position_dodge(0.7), size=0.4, width = 0.45) +
  xlab("Treatment") +
  ylab("Weight per seed (mg)") +
  facet_wrap(~ Species, ncol =3, scales = "free_y", labeller=label_parsed) + # need the labeller function to parse the italic text
  scale_fill_manual(values = vir_lite(viridis(5), ds=0.6, dv=0.6), name = "Species") + 
  scale_x_discrete(labels = c("Control", "Heat", "Heat+\nWater", "Water")) +
  guides(fill = guide_legend(label.theme = element_text(face = "italic", size = 9), ncol=3)) +
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.05))) +
  theme_bw() +
  theme(legend.title=element_text(size=10),
        legend.text = element_text(size=9.5),
        axis.text = element_text(colour = "black"),
        panel.spacing.y = unit(0.2, "lines"),
        panel.spacing.x = unit(0.4, "lines"),
        panel.grid.major.x = element_blank(),
        strip.text.x = element_text(size = 10, margin = margin(0,0,0,0, "mm")),
        axis.text.x=element_text(size=9.5, margin=margin(t=1, r=0, b=0, l=0)),
        axis.title.x=element_text(size=10, margin=margin(t=-2, r=0, b=0, l=0)),
        axis.title.y=element_text(size=10, margin=margin(t=0, r=1, b=0, l=0)),
        legend.margin = margin(t=15, r=0, b=0, l=0),
        plot.margin=unit(c(1,1,1,1),"mm"))

plot_panel = plot_grid(sn_plot,
                       lemon::reposition_legend(sw_plot, 'center', panel=c('panel-3-2','panel-3-3'), plot=TRUE),
                       labels = c("(a)", "(b)"),
                       label_size = 11, #scale = 1, label_colour = "blue",
                       #label_x = 0, label_y = 0.11,
                       #vjust = -0.005,
                       align = "v", axis = 'l',
                       ncol = 1,
                       rel_heights = c(1, 1.39))


ppi = 300
png("Figures/Seeds_cowplot_7-8.png", width=7*ppi, height=8*ppi, res=ppi)
#pdf("Figures/Seeds_cowplot_7-8.pdf", width=7, height=8)
plot_grid(plot_panel)
dev.off()



#
#####