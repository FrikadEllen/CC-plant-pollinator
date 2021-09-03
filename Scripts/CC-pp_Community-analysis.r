.libPaths(c("E:\\Rlib", .libPaths()))

library(vegan)
library(pairwiseAdonis)


#### Data

raw_int=read.csv("Data/Raw_InteractionData_BothYears.csv", header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
raw_int$Record_ID = 1:length(raw_int$Specimen_ID)
raw_int$Plot2 = ifelse(raw_int$Year == 2014, raw_int$Plot, raw_int$Plot + 24)

raw_int$Round = factor(raw_int$Round)
raw_int$Specimen_ID = factor(raw_int$Specimen_ID)

str(raw_int)
# 'data.frame':	3882 obs. of  9 variables:
# $ Year       : int  2014 2014 2014 2014 2014 2014 2014 2014 2014 2014 ...
# $ Plot       : int  2 2 16 1 11 23 23 9 9 15 ...
# $ Round      : Factor w/ 9 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Plant      : chr  "Veronica_persica" "Lamium_purpureum" "Stellaria_media" "Stellaria_media" ...
# $ Specimen_ID: Factor w/ 104 levels "Aglais_io","Aglais_urticae",..: 97 17 47 8 47 103 98 8 8 98 ...
# $ Basic_Type : chr  "Hoverfly" "Bumblebee" "Hoverfly" "Fly" ...
# $ Treatment  : chr  "Heat+Water" "Heat+Water" "Control" "Control" ...
# $ Record_ID  : int  1 2 3 4 5 6 7 8 9 10 ...
# $ Plot2      : num  2 2 16 1 11 23 23 9 9 15 ...
head(raw_int)
#   Year Plot Round            Plant      Specimen_ID Basic_Type  Treatment Record_ID Plot2
# 1 2014    2     1 Veronica_persica Sphaerophoria_sp   Hoverfly Heat+Water         1     2
# 2 2014    2     1 Lamium_purpureum Bombus_pascuorum  Bumblebee Heat+Water         2     2
# 3 2014   16     1  Stellaria_media  Eristalis_tenax   Hoverfly    Control         3    16
# 4 2014    1     1  Stellaria_media  Anthomyiidae_sp        Fly    Control         4     1
# 5 2014   11     1 Veronica_persica  Eristalis_tenax   Hoverfly       Heat         5    11
# 6 2014   23     1 Veronica_persica    Tachinidae_sp        Fly    Control         6    23
tail(raw_int)
#      Year Plot Round                 Plant            Specimen_ID Basic_Type Treatment Record_ID Plot2
# 3877 2015   24     7 Chrysanthemum_segetum        Syritta_pipiens   Hoverfly     Water      3877    48
# 3878 2015   24     7 Chrysanthemum_segetum  Platycheirus_occultus   Hoverfly     Water      3878    48
# 3879 2015   24     7 Chrysanthemum_segetum       Lucilia_silvarum        Fly     Water      3879    48
# 3880 2015   24     7 Chrysanthemum_segetum       Eupeodes_luniger   Hoverfly     Water      3880    48
# 3881 2015   24     7 Chrysanthemum_segetum Platycheirus_albimanus   Hoverfly     Water      3881    48
# 3882 2015   24     7 Chrysanthemum_segetum        Eristalis_tenax   Hoverfly     Water      3882    48


insects_mat = as.data.frame(tapply(raw_int$Record_ID, list(raw_int$Plot2, raw_int$Specimen_ID), length))
insects_mat[is.na(insects_mat)] = 0

flowers_mat = read.csv("Data/Raw_FlowerMatrix_BothYears.csv", header=TRUE, sep=",", as.is=TRUE, na.strings="(null)", row.names = 1)

years2=read.table("Data/SimilarityIndex/TreatmentYearData.txt", h=T, dec=".")
years2$Treatment = as.factor(years2$Treatment)
years2$Year = as.factor(years2$Year)





# 1. Flowering plants ~ year and treatment 
# 2. Extrapolated Insects ~ year and treatment



########################### 1. Flowering plants ~ year and treatment ####

adonis2(flowers_mat~years2$Treatment*years2$Year, method="bray", binary=TRUE, permutations=999, by = "terms")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = flowers_mat ~ years2$Treatment * years2$Year, permutations = 999, method = "bray", by = "terms", binary = TRUE)
#                              Df SumOfSqs      R2       F Pr(>F)    
# years2$Treatment              3  0.16470 0.07500  1.7348  0.059 .  
# years2$Year                   1  0.65682 0.29911 20.7546  0.001 ***
# years2$Treatment:years2$Year  3  0.10855 0.04943  1.1434  0.365    
# Residual                     40  1.26588 0.57646                   
# Total                        47  2.19595 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis2(flowers_mat~years2$Treatment+years2$Year, method="bray", binary=TRUE, permutations=999, by = "margin")
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = flowers_mat ~ years2$Treatment + years2$Year, permutations = 999, method = "bray", by = "margin", binary = TRUE)
#                  Df SumOfSqs      R2       F Pr(>F)    
# years2$Treatment  3  0.16470 0.07500  1.7176  0.074 .  
# years2$Year       1  0.65682 0.29911 20.5491  0.001 ***
# Residual         43  1.37443 0.62589                   
# Total            47  2.19595 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##### validation ####

fdist1=vegdist(flowers_mat, method="bray", binary=TRUE)
fsim=anosim(fdist1, years2$Treatment)
summary(fsim)
# Call:
#   anosim(x = fdist1, grouping = years2$Treatment) 
# Dissimilarity: binary bray 
# 
# ANOSIM statistic R: 0.01017 
# Significance: 0.354 
# 
# Permutation: free
# Number of permutations: 999
# 
# Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99% 
#   0.0427 0.0578 0.0746 0.0923 
# 
# Dissimilarity ranks between and within classes:
#   0%     25%    50%    75%   100%   N
# Between     1.0 293.000 581.50 857.00 1127.0 864
# Control     6.0 325.000 573.75 818.75 1006.0  66
# Heat       13.0 513.500 774.50 993.00 1127.0  66
# Heat+Water 38.5 364.000 589.25 880.50 1057.0  66
# Water       4.5 180.375 364.00 654.50 1068.5  66

plot(fsim)


# Betadisperser: test of equality of variances
test1=with(years2, betadisper(fdist1, years2$Treatment))
test1
# Homogeneity of multivariate dispersions 
# Call: betadisper(d = fdist1, group = years2$Treatment) 
# No. of Positive Eigenvalues: 19
# No. of Negative Eigenvalues: 28 

# Average distance to median:
#   Control       Heat Heat+Water      Water 
# 0.1929     0.2313     0.1971     0.1595 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 47 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 0.8601 0.4227 0.3388 0.2881 0.2052 0.1862 0.1505 0.1240 

plot(test1)
boxplot(test1)

anova(test1)
# Analysis of Variance Table 
# Response: Distances
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     3 0.030964 0.0103214    2.36 0.08438 .
# Residuals 44 0.192434 0.0043735                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutest(test1)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq    F N.Perm Pr(>F)  
# Groups     3 0.030964 0.0103214 2.36    999  0.089 .
# Residuals 44 0.192434 0.0043735                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


TukeyHSD(test1)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                            diff         lwr         upr     p adj
# Heat-Control        0.038368317 -0.03371774 0.110454370 0.4932377
# Heat+Water-Control  0.004170752 -0.06791530 0.076256806 0.9986662
# Water-Control      -0.033345908 -0.10543196 0.038740146 0.6082802
# Heat+Water-Heat    -0.034197564 -0.10628362 0.037888489 0.5886217
# Water-Heat         -0.071714225 -0.14380028 0.000371829 0.0516439
# Water-Heat+Water   -0.037516660 -0.10960271 0.034569393 0.5124594



#
########################### 2. Extrapolated Insects ~ year and treatment ####

adonis2(insects_mat~years2$Treatment*years2$Year, method="chao", permutations=999, by = "terms")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = insects_mat ~ years2$Treatment * years2$Year, permutations = 999, method = "chao", by = "terms")
#                              Df SumOfSqs      R2       F Pr(>F)    
# years2$Treatment              3  0.19218 0.12172  3.9034  0.008 ** 
# years2$Year                   1  0.70332 0.44548 42.8563  0.001 ***
# years2$Treatment:years2$Year  3  0.02684 0.01700  0.5452  0.735    
# Residual                     40  0.65645 0.41579                   
# Total                        47  1.57879 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


adonis2(insects_mat~years2$Treatment+years2$Year, method="chao", permutations=999, by = "margin")
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = insects_mat ~ years2$Treatment + years2$Year, permutations = 999, method = "chao", by = "margin")
#                  Df SumOfSqs      R2       F Pr(>F)    
# years2$Treatment  3  0.19218 0.12172  4.0313  0.004 ** 
# years2$Year       1  0.70332 0.44548 44.2606  0.001 ***
# Residual         43  0.68329 0.43279                   
# Total            47  1.57879 1.00000                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


pairwise.adonis2(insects_mat ~ Treatment + Year, method="chao", permutations=999, by = "margin", data = years2)
# $parent_call
# [1] "insects_mat ~ Treatment + Year , strata = Null"
# 
# $`Control_vs_Heat+Water`
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.10488 0.10488  4.1014 0.10614  0.021 *  
# Year       1   0.34625 0.34625 13.5403 0.35041  0.001 ***
# Residuals 21   0.53700 0.02557         0.54345           
# Total     23   0.98813                 1.00000           
# ---
# 
# $Control_vs_Water
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.01545 0.01545   1.281 0.02335  0.334    
# Year       1   0.39296 0.39296  32.569 0.59379  0.001 ***
# Residuals 21   0.25337 0.01207         0.38286           
# Total     23   0.66179                 1.00000           
# ---
# 
# $Control_vs_Heat
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.03784 0.03784   3.215 0.05084  0.067 .  
# Year       1   0.45934 0.45934  39.029 0.61712  0.001 ***
# Residuals 21   0.24716 0.01177         0.33205           
# Total     23   0.74434                 1.00000           
# ---
# 
# $`Heat+Water_vs_Water`
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.14216 0.142160  7.9091 0.17289  0.001 ***
# Year       1   0.30265 0.302646 16.8376 0.36806  0.001 ***
# Residuals 21   0.37746 0.017974         0.45905           
# Total     23   0.82227                  1.00000           
# ---
# 
# $`Heat+Water_vs_Heat`
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.02839 0.02839  1.4303 0.03694  0.306    
# Year       1   0.32342 0.32342 16.2930 0.42075  0.001 ***
# Residuals 21   0.41686 0.01985         0.54231           
# Total     23   0.76867                 1.00000           
# ---
# 
# $Water_vs_Heat
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Treatment  1   0.05563 0.055631   6.111 0.09952  0.015 *  
# Year       1   0.31219 0.312188  34.293 0.55848  0.001 ***
# Residuals 21   0.19118 0.009104         0.34200           
# Total     23   0.55899                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1       


pairwise.adonis2(insects_mat ~ Treatment, method="chao", permutations=999, by = "margin", data = years2)
# $parent_call
# [1] "insects_mat ~ Treatment , strata = Null"
# 
# $`Control_vs_Heat+Water`
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Treatment  1   0.10488 0.104879  2.6123 0.10614  0.078 .
# Residuals 22   0.88325 0.040148         0.89386         
# Total     23   0.98813                  1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $Control_vs_Water
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Treatment  1   0.01545 0.015450  0.5259 0.02335  0.612
# Residuals 22   0.64633 0.029379         0.97665       
# Total     23   0.66179                  1.00000       
# 
# $Control_vs_Heat
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Treatment  1   0.03784 0.037840  1.1783 0.05084  0.358
# Residuals 22   0.70650 0.032114         0.94916       
# Total     23   0.74434                  1.00000       
# 
# $`Heat+Water_vs_Water`
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# Treatment  1   0.14216 0.142160  4.5986 0.17289  0.014 *
# Residuals 22   0.68011 0.030914         0.82711         
# Total     23   0.82227                  1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $`Heat+Water_vs_Heat`
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Treatment  1   0.02839 0.028392 0.84377 0.03694   0.49
# Residuals 22   0.74028 0.033649         0.96306       
# Total     23   0.76867                  1.00000       
# 
# $Water_vs_Heat
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Treatment  1   0.05563 0.055631  2.4314 0.09952  0.146
# Residuals 22   0.50336 0.022880         0.90048       
# Total     23   0.55899                  1.00000              
  


##### Validation ####

fdist1=vegdist(insects_mat, method="chao")

# Betadisperser: test of equality of variances
test1=with(years2, betadisper(fdist1, years2$Treatment))
# Warning message:
#   In betadisper(fdist1, years2$Treatment) :
#   some squared distances are negative and changed to zero
test1
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = fdist1, group = years2$Treatment)
# 
# No. of Positive Eigenvalues: 23
# No. of Negative Eigenvalues: 24
# 
# Average distance to median:
#   Control       Heat Heat+Water      Water 
# 0.1589     0.1380     0.1684     0.1103 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 47 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 1.0142 0.4213 0.3082 0.2635 0.2185 0.2045 0.1790 0.1611 

plot(test1)
boxplot(test1)

anova(test1)
# Analysis of Variance Table
# 
# Response: Distances
#           Df  Sum Sq   Mean Sq F value Pr(>F)
# Groups     3 0.02386 0.0079541  0.8053 0.4977
# Residuals 44 0.43462 0.0098778                    

permutest(test1)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     3 0.02386 0.0079541 0.8053    999  0.498
# Residuals 44 0.43462 0.0098778   

TukeyHSD(test1)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
# diff         lwr         upr     p adj
# Heat-Control       -0.020829843 -0.12916419 0.08750451 0.9554667
# Heat+Water-Control  0.009581697 -0.09875265 0.11791605 0.9952958
# Water-Control      -0.048545499 -0.15687985 0.05978885 0.6323061
# Heat+Water-Heat     0.030411540 -0.07792281 0.13874589 0.8763760
# Water-Heat         -0.027715656 -0.13605001 0.08061869 0.9029861
# Water-Heat+Water   -0.058127196 -0.16646155 0.05020715 0.4862970


test2=with(years2, betadisper(fdist1, years2$Year))
test2
# Homogeneity of multivariate dispersions
# 
# Call: betadisper(d = fdist1, group = years2$Year)
# 
# No. of Positive Eigenvalues: 23
# No. of Negative Eigenvalues: 24
# 
# Average distance to median:
#   1       2 
# 0.13592 0.09405 
# 
# Eigenvalues for PCoA axes:
#   (Showing 8 of 47 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 1.0142 0.4213 0.3082 0.2635 0.2185 0.2045 0.1790 0.1611 

plot(test2)
boxplot(test2)

anova(test2)
# Analysis of Variance Table
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq F value  Pr(>F)  
# Groups     1 0.021045 0.0210446   3.837 0.05621 .
# Residuals 46 0.252291 0.0054846                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutest(test2)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)  
# Groups     1 0.021045 0.0210446 3.837    999  0.056 .
# Residuals 46 0.252291 0.0054846                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#
#####