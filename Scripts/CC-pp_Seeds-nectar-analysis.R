.libPaths(c("E:\\Rlib", .libPaths()))


library(lme4)
library(emmeans)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(plotrix)


######################### Functions


######################### Datasets


# 1. 2014 Cornflower Seed Number
# 2. 2014 Cornflower Average Seed Weight
# 3. 2014 Corn Marigold Seed Number
# 4. 2014 Corn Marigold Average Seed Weight
# 5. 2015 Cornflower Seed Number
# 6. 2015 Cornflower Average Seed Weight
# 7. 2015 Corn Marigold Seed Number
# 8. 2015 Corn Marigold Average Seed Weight
# 9. 2015 Red Dead Nettle Average Seed Weight 
# 10. 2015 Field Speedwell Average Seed Number 
# 11. 2015 Field Speedwell Average Seed Weight 
# 12. 2015 Common Chickweed Average Seed Number 
# 13. 2015 Common Chickweed Average Seed Weight
# 14. 2015 Cornflower Nectar Volume
# 15. 2015 Red Dead Nettle Nectar Volume
# 16. 2015 Speedwell Nectar Volume
# 17. Corn Marigold Flower Diameter


# load the data file:
Seeds = read.csv("Data/AllData_MultiPlot_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(Seeds)
# change plot to factor:
Seeds$Plot = as.factor(Seeds$Plot)
# change Treatment to factor:
Seeds$Treatment = as.factor(Seeds$Treatment)
# some of the columns have been classed as character, need to set all relevant ones to numeric:
Seeds[, c(3:4, 6:8, 10:22, 24:27)] = lapply(Seeds[, c(3:4, 6:8, 10:22, 24:27)], as.numeric)

# change the seed weight columns from grams to mg (so the values are easier to interpret):
Seeds$CF14_AWeight = Seeds$CF14_AWeight *1000
Seeds$CM14_AWeight = Seeds$CM14_AWeight*1000
Seeds$CF15_AWeight = Seeds$CF15_AWeight *1000
Seeds$CM15_AWeight = Seeds$CM15_AWeight*1000
Seeds$RDN15_AWeight = Seeds$RDN15_AWeight*1000
Seeds$SW15_AWeight = Seeds$SW15_AWeight*1000
Seeds$CW15_AWeight = Seeds$CW15_AWeight*1000


# N.B. the "XXX_Day" columns are the day of the year of the associated date. e.g.
# CM14_Date	  CM14_Day
# 31/07/2014	212
# so this date is the 212th day of the year.



#
########################## 1. 2014 Cornflower Seed Number ####

# quick look at the distribution:
plot(density(Seeds$CF14_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glmm1=glmer(CF14_Seeds~Treatment +(1|Plot) + (1|CF14_Day), data=Seeds, family=poisson, na.action = "na.omit")
isSingular(glmm1)
# TRUE - the variance of Plot is 0
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment + CF14_Day , abline = 0)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.058562, p-value = 0.3828
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.93312, p-value = 0.168

# negative binomial:
glmm2=glmer.nb(CF14_Seeds~Treatment +(1|Plot) + (1|CF14_Day), data=Seeds, na.action = "na.omit")
# Warning message:
#   In theta.ml(Y, mu, weights = object@resp$weights, limit = limit,  :
#                 iteration limit reached
# this model struggled
isSingular(glmm2)
# [1] FALSE - but the variance of Plot is really small (1.877e-10)
plot(glmm2)
plot(glmm2, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
plot(glmm2, resid(., scaled=TRUE) ~ fitted(.) | Treatment + CF14_Day , abline = 0)
fitnb <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnb)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.051042, p-value = 0.5593
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.93694, p-value = 0.176

# Poisson looks the best, but I will drop Plot:
glmm1a=glmer(CF14_Seeds~Treatment + (1|CF14_Day), data=Seeds, family=poisson, na.action = "na.omit")
plot(glmm1a)
fitpois <- simulateResiduals(fittedModel = glmm1a, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.049666, p-value = 0.5946
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.9334, p-value = 0.184


summary(glmm1a)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: CF14_Seeds ~ Treatment + (1 | CF14_Day)
# Data: Seeds
# 
# AIC      BIC   logLik deviance df.resid 
# 1414.2   1431.6   -702.1   1404.2      235 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.98448 -0.60744  0.00249  0.60521  2.89431 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# CF14_Day (Intercept) 0.001784 0.04224 
# Number of obs: 240, groups:  CF14_Day, 2
# 
# Fixed effects:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          3.23123    0.03937  82.071  < 2e-16 ***
# TreatmentHeat       -0.17733    0.03799  -4.668 3.04e-06 ***
# TreatmentHeat+Water -0.15942    0.03781  -4.217 2.48e-05 ***
# TreatmentWater      -0.05337    0.03677  -1.452    0.147    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.440              
# TrtmntHt+Wt -0.442  0.458       
# TreatmntWtr -0.454  0.471  0.473

drop1(glmm1a, test="Chisq")
# Single term deletions
# 
# Model:
#   CF14_Seeds ~ Treatment + (1 | CF14_Day)
#           Df    AIC    LRT   Pr(Chi)    
# <none>       1414.2                     
# Treatment  3 1438.3 30.125 1.299e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glmm1a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment    lsmean         SE df asymp.LCL asymp.UCL
# Control    3.231232 0.03937129 NA  3.154065  3.308398
# Heat       3.053897 0.04096015 NA  2.973617  3.134178
# Heat+Water 3.071807 0.04078960 NA  2.991861  3.151753
# Water      3.177858 0.03982665 NA  3.099799  3.255917
# 
# Results are given on the log scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate         SE df    z.ratio p.value
# Control - Heat        0.17733450 0.03799028 NA  4.6678916  <.0001
# Control - Heat+Water  0.15942466 0.03780634 NA  4.2168769  0.0001
# Control - Water       0.05337378 0.03676533 NA  1.4517425  0.4670
# Heat - Heat+Water    -0.01790984 0.03945826 NA -0.4538933  0.9689
# Heat - Water         -0.12396072 0.03846199 NA -3.2229409  0.0070
# Heat+Water - Water   -0.10605088 0.03828031 NA -2.7703764  0.0286
# 
# Results are given on the log scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CF14_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Seeds per Cornflower Seed Head",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CF14_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CF14_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 25.33333   21.21667   21.60000   24.01667 
prse
#   Control       Heat Heat+Water      Water 
# 0.5859305  0.5518800  0.5665836  0.6437407
prbarplot=barplot(prmean, ylim=c(0,30), xlab="Treatment",ylab="Mean Seeds per Cornflower Seed Head",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 2. 2014 Cornflower Average Seed Weight ####

# quick look at the distribution:
plot(density(Seeds$CF14_AWeight, na.rm=TRUE))

mm1 = lmer(CF14_AWeight ~ Treatment +(1|Plot) + (1|CF14_Day), data=Seeds, REML=F, na.action = "na.omit") 
# warning message about boundary:
isSingular(mm1)
# [1] FALSE
shapiro.test(residuals(mm1))
# data:  residuals(mm1)
# W = 0.99229, p-value = 0.2434
qqnorm(residuals(mm1))
qqline(residuals(mm1))
plot(mm1)

# REML = TRUE:
mm1a = lmer(CF14_AWeight ~ Treatment +(1|Plot) + (1|CF14_Day), data=Seeds, REML=T, na.action = "na.omit")
isSingular(mm1a)
# [1] FALSE
shapiro.test(residuals(mm1a))
# data:  residuals(mm1a)
# W = 0.99274, p-value = 0.2896
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm1a, n = 250)
plot(fitg, asFactor = TRUE)

# this looks pretty good, but there is some right skew so it's worth trying Gamma

# Gamma:
glmm1=glmer(CF14_AWeight ~ Treatment +(1|Plot) + (1|CF14_Day), data=Seeds, family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, n = 250, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.6106, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.13552, p-value < 2.2e-16

# gaussian is a better fit

summary(mm1a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CF14_AWeight ~ Treatment + (1 | Plot) + (1 | CF14_Day)
# Data: Seeds
# 
# REML criterion at convergence: 623.7
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.38300 -0.63962 -0.06072  0.61359  2.98172 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.009378 0.09684 
# CF14_Day (Intercept) 0.332105 0.57629 
# Residual             0.747148 0.86438 
# Number of obs: 240, groups:  Plot, 24; CF14_Day, 2
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)          3.49897    0.42434   8.246
# TreatmentHeat        0.24285    0.16742   1.451
# TreatmentHeat+Water  0.04555    0.16742   0.272
# TreatmentWater      -0.33001    0.16742  -1.971
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.197              
# TrtmntHt+Wt -0.197  0.500       
# TreatmntWtr -0.197  0.500  0.500

drop1(mm1a, test="Chisq")
# Single term deletions
# 
# Model:
#   CF14_AWeight ~ Treatment + (1 | Plot) + (1 | CF14_Day)
#           npar    AIC    LRT Pr(Chi)  
# <none>         631.25                 
# Treatment    3 636.14 10.895 0.01231 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


lsmeans(mm1a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean    SE   df lower.CL upper.CL
# Control      3.50 0.424 1.13   -0.607     7.61
# Heat         3.74 0.424 1.13   -0.365     7.85
# Heat+Water   3.54 0.424 1.13   -0.562     7.65
# Water        3.17 0.424 1.13   -0.937     7.28
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE df t.ratio p.value
# Control - Heat          -0.2429 0.167 20 -1.451  0.4843 
# Control - (Heat+Water)  -0.0455 0.167 20 -0.272  0.9927 
# Control - Water          0.3300 0.167 20  1.971  0.2318 
# Heat - (Heat+Water)      0.1973 0.167 20  1.178  0.6468 
# Heat - Water             0.5729 0.167 20  3.422  0.0132 
# (Heat+Water) - Water     0.3756 0.167 20  2.243  0.1458 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CF14_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Cornflower Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CF14_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CF14_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 3.498969   3.741822   3.544514   3.168955 
prse
# Control       Heat Heat+Water      Water 
# 0.1370526  0.1144944  0.1323163  0.1104420 
prbarplot=barplot(prmean, 
                  ylim=c(0,0.004), 
                  xlab="Treatment",ylab="Mean Cornflower Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 3. 2014 Corn Marigold Seed Number ####

# quick look at the distribution:
plot(density(Seeds$CM14_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glmm1=glmer(CM14_Seeds ~ Treatment +(1|Plot) + (1|CM14_Day), data=Seeds[-which(is.na(Seeds$CM14_Seeds)),], family=poisson, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.077293, p-value = 0.3818
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.2405, p-value = 0.112

# negative binomial:
glmm2=glmer.nb(CM14_Seeds ~ Treatment +(1|Plot) + (1|CM14_Day), data=Seeds[-which(is.na(Seeds$CM14_Seeds)),], na.action = "na.omit")
isSingular(glmm2)
# TRUE - the variance of Day is almost 0 (5.094e-11)
plot(glmm2)
plot(glmm2, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitnb <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnb)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.086202, p-value = 0.2567
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.89011, p-value = 0.232

# Poisson looks the best fit

summary(glmm1)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: poisson  ( log )
# Formula: CM14_Seeds ~ Treatment + (1 | Plot) + (1 | CM14_Day)
# Data: Seeds
# 
# AIC      BIC   logLik deviance df.resid 
# 2066.3   2083.9  -1027.1   2054.3      132 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -7.5136 -1.5403  0.1438  1.7408  7.9109 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# Plot     (Intercept) 2.462e-02 0.15689 
# CM14_Day (Intercept) 2.218e-05 0.00471 
# Number of obs: 138, groups:  Plot, 23; CM14_Day, 2
# 
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          5.6307390  0.0649106  86.746   <2e-16 ***
# TreatmentHeat       -0.2275837  0.0918472  -2.478   0.0132 *  
# TreatmentHeat+Water -0.2170725  0.0963053  -2.254   0.0242 *  
# TreatmentWater       0.0001551  0.0916764   0.002   0.9987    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.705              
# TrtmntHt+Wt -0.672  0.475       
# TreatmntWtr -0.706  0.499  0.476

drop1(glmm1, test="Chisq")
# Single term deletions
# 
# Model:
#   CM14_Seeds ~ Treatment + (1 | Plot) + (1 | CM14_Day)
#           Df    AIC    LRT Pr(Chi)  
# <none>       2066.3                 
# Treatment  3 2069.5 9.1858 0.02692 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glmm1, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment    lsmean         SE df asymp.LCL asymp.UCL
# Control    5.630739 0.06491060 NA  5.503517  5.757961
# Heat       5.403155 0.06515147 NA  5.275461  5.530850
# Heat+Water 5.413666 0.07129873 NA  5.273924  5.553409
# Water      5.630894 0.06491037 NA  5.503672  5.758116
# 
# Results are given on the log scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                  estimate         SE df      z.ratio p.value
# Control - Heat        0.2275836521 0.09184718 NA  2.477851307  0.0634
# Control - Heat+Water  0.2170724878 0.09630532 NA  2.254003068  0.1090
# Control - Water      -0.0001550718 0.09167639 NA -0.001691513  1.0000
# Heat - Heat+Water    -0.0105111644 0.09646694 NA -0.108961311  0.9995
# Heat - Water         -0.2277387239 0.09184703 NA -2.479543649  0.0631
# Heat+Water - Water   -0.2172275595 0.09630515 NA -2.255617196  0.1086
# 
# Results are given on the log scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CM14_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Mean Corn Marigold Seed Number",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CM14_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CM14_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 279.0833   230.6389   225.7333   279.0833 
prse
# Control       Heat Heat+Water      Water 
# 7.115975  12.294479   8.181270   7.285164 
prbarplot=barplot(prmean, 
                  ylim=c(0,300), 
                  xlab="Treatment",ylab="Mean Corn Marigold Seed Number",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)




########################## 4. 2014 Corn Marigold Average Seed Weight ####

# quick look at the distribution:
plot(density(Seeds$CM14_AWeight, na.rm=TRUE))

mm1 = lmer(CM14_AWeight ~ Treatment +(1|Plot) + (1|CM14_Day), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# data:  residuals(mm1)
# W = 0.99129, p-value = 0.5522
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)
plot(mm1)

# REML = TRUE:
mm1a = lmer(CM14_AWeight ~ Treatment +(1|Plot) + (1|CM14_Day), data=Seeds[-which(is.na(Seeds$CM14_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# data:  residuals(mm1a)
# W = 0.99257, p-value = 0.6876
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

fitg <- simulateResiduals(fittedModel = mm1a, n = 250)
plot(fitg, asFactor = TRUE)

# This model looks OK, but worth trying gamma

# Gamma:
glmm1=glmer(CM14_AWeight ~ Treatment +(1|Plot) + (1|CM14_Day), data=Seeds[-which(is.na(Seeds$CM14_AWeight)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.64131, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.13089, p-value < 2.2e-16

# gaussian is the best model

summary(mm1a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CM14_AWeight ~ Treatment + (1 | Plot) + (1 | CM14_Day)
# Data: Seeds[-which(is.na(Seeds$CM14_AWeight)), ]
# 
# REML criterion at convergence: 129.4
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.51224 -0.67616 -0.09983  0.73982  3.00063 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.030907 0.17580 
# CM14_Day (Intercept) 0.003046 0.05519 
# Residual             0.120342 0.34690 
# Number of obs: 138, groups:  Plot, 23; CM14_Day, 2
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)           1.6681     0.1001  16.667
# TreatmentHeat        -0.3487     0.1303  -2.675
# TreatmentHeat+Water  -0.5257     0.1367  -3.846
# TreatmentWater       -0.0956     0.1303  -0.733
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.651              
# TrtmntHt+Wt -0.621  0.477       
# TreatmntWtr -0.651  0.500  0.477

drop1(mm1a, test="Chisq")
# Single term deletions 
# Model:
#   CM14_AWeight ~ Treatment + (1 | Plot) + (1 | CM14_Day)
#           npar    AIC    LRT  Pr(Chi)   
# <none>         131.91                   
# Treatment    3 141.45 15.543 0.001407 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(mm1a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean    SE   df lower.CL upper.CL
# Control      1.67 0.100 10.5    1.447     1.89
# Heat         1.32 0.100 10.5    1.098     1.54
# Heat+Water   1.14 0.108 12.2    0.907     1.38
# Water        1.57 0.100 10.5    1.351     1.79
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE df t.ratio p.value
# Control - Heat           0.3487 0.130 19  2.675  0.0656 
# Control - (Heat+Water)   0.5257 0.137 19  3.846  0.0055 
# Control - Water          0.0956 0.130 19  0.733  0.8824 
# Heat - (Heat+Water)      0.1770 0.137 19  1.295  0.5771 
# Heat - Water            -0.2531 0.130 19 -1.942  0.2445 
# (Heat+Water) - Water    -0.4301 0.137 19 -3.146  0.0251 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CM14_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CM14_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CM14_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 1.668132   1.319430   1.142424   1.572530 
prse
# Control       Heat Heat+Water      Water 
# 0.07718809 0.06372171 0.06181955 0.05571544
prbarplot=barplot(prmean, 
                  ylim=c(0,0.002), 
                  xlab="Treatment",ylab="Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 5. 2015 Cornflower Seed Number ####

# quick look at the distribution:
plot(density(Seeds$CF15_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# only one collection date for cornflower in 2015, so only 1 random effect (plot)

# poisson:
glmm1=glmer(CF15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Seeds)),], family=poisson, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.13184, p-value = 0.2947
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.2928, p-value = 0.024

# poisson is overdispersed, try negative binomial:
glmm2=glmer.nb(CF15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Seeds)),], na.action = "na.omit")
isSingular(glmm2)
# [1] TRUE - variance of plot is 0
plot(glmm2)
plot(glmm2, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitnb <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnb)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.12951, p-value = 0.3148
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.87527, p-value = 0.288

# neg bin is better, but is singular, drop Plot from the model:
glm1 = glm.nb(CF15_Seeds ~ Treatment, data=Seeds[-which(is.na(Seeds$CF15_Seeds)),], na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
fitnb2 <- simulateResiduals(fittedModel = glm1, n = 250)
plot(fitnb2)
testUniformity(fitnb2)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.13741, p-value = 0.2501
testDispersion(fitnb2)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.88784, p-value = 0.312

# check that the random effect (Plot) really isn't useful:
pchisq(2*(logLik(glmm2)-logLik(glm1)),
       df=1,lower.tail=FALSE)/2
# 'log Lik.' 0.5 (df=6)

summary(glm1)
#   glm.nb(formula = CF15_Seeds ~ Treatment, data = Seeds, na.action = "na.omit", 
#          init.theta = 10.77813001, link = log) 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.3589  -0.4722   0.1063   0.5809   1.4375  
# 
# Coefficients:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          3.12613    0.09880  31.640   <2e-16 ***
# TreatmentHeat       -0.17048    0.14180  -1.202    0.229    
# TreatmentHeat+Water -0.25445    0.14043  -1.812    0.070 .  
# TreatmentWater       0.04145    0.14492   0.286    0.775    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(10.7781) family taken to be 1)
# 
# Null deviance: 71.297  on 54  degrees of freedom
# Residual deviance: 65.614  on 51  degrees of freedom
# (185 observations deleted due to missingness)
# AIC: 395.02
# 
# Number of Fisher Scoring iterations: 1 
# 
# Theta:  10.78 
# Std. Err.:  3.59  
# 2 x log-likelihood:  -385.021 

anova(glm1, test="Chisq")
# Analysis of Deviance Table 
# Model: Negative Binomial(10.7781), link: log 
# Response: CF15_Seeds 
# Terms added sequentially (first to last) 
# 
# Df Deviance Resid. Df Resid. Dev Pr(>Chi)
# NULL                         54     71.297         
# Treatment  3   5.6826        51     65.614   0.1281

drop1(glm1, test="Chisq")
# Single term deletions 
# Model:
#   CF15_Seeds ~ Treatment
#           Df Deviance    AIC    LRT Pr(>Chi)
# <none>         65.614 393.02                
# Treatment  3   71.297 392.70 5.6826   0.1281

lsmeans(glm1, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment    lsmean         SE df asymp.LCL asymp.UCL
# Control    3.126134 0.09880270 NA  2.932484  3.319784
# Heat       2.955654 0.10170865 NA  2.756309  3.154999
# Heat+Water 2.871680 0.09979454 NA  2.676086  3.067273
# Water      3.167583 0.10602112 NA  2.959785  3.375380
# 
# Results are given on the log scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate        SE df    z.ratio p.value
# Control - Heat        0.17047972 0.1417978 NA  1.2022732  0.6254
# Control - Heat+Water  0.25445415 0.1404312 NA  1.8119488  0.2676
# Control - Water      -0.04144876 0.1449222 NA -0.2860069  0.9919
# Heat - Heat+Water     0.08397443 0.1424907 NA  0.5893327  0.9353
# Heat - Water         -0.21192848 0.1469188 NA -1.4424874  0.4727
# Heat+Water - Water   -0.29590291 0.1456002 NA -2.0322969  0.1761
# 
# Results are given on the log scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CF15_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Mean Cornflower Seed Number",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CF15_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CF15_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 22.78571   19.21429   17.66667   23.75000 
prse
# Control       Heat Heat+Water      Water 
# 1.8282838  0.9089214  2.0087112  2.6489706 
prbarplot=barplot(prmean, ylim=c(0,30), xlab="Treatment",ylab="Mean Cornflower Seed Number",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 6. 2015 Cornflower Average Seed Weight ####

# quick look at the distribution:
plot(density(Seeds$CF15_AWeight, na.rm=TRUE))

# only one collection date for cornflower in 2015, so only 1 random effect (plot)
mm1 = lmer(CF15_AWeight ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# data:  residuals(mm1)
# W = 0.96351, p-value = 0.09344
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)
plot(mm1)

# REML = TRUE:
mm1a = lmer(CF15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# data:  residuals(mm1a)
# W = 0.96617, p-value = 0.1235
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

# The gaussian model looks OK, but worth trying gamma

glmm1=glmer(CF15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_AWeight)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, n=250, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.65529, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.15714, p-value < 2.2e-16

# gaussian is definitely the best fit

summary(mm1a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CF15_AWeight ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$CF15_AWeight)), ]
# 
# REML criterion at convergence: 147.9
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.60980 -0.52898 -0.05705  0.47343  2.16781 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.3619   0.6015  
# Residual             0.6506   0.8066  
# Number of obs: 55, groups:  Plot, 20
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)           3.4218     0.3318  10.312
# TreatmentHeat         0.2195     0.4792   0.458
# TreatmentHeat+Water  -0.3504     0.4752  -0.737
# TreatmentWater       -0.7835     0.5048  -1.552
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.692              
# TrtmntHt+Wt -0.698  0.483       
# TreatmntWtr -0.657  0.455  0.459

drop1(mm1a, test="Chisq")
# Single term deletions 
# Model:
#   CF15_AWeight ~ Treatment + (1 | Plot)
#           npar    AIC   LRT Pr(Chi)
# <none>         158.32              
# Treatment    3 157.21 4.892  0.1799

lsmeans(mm1a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean    SE   df lower.CL upper.CL
# Control      3.42 0.333 17.9     2.72     4.12
# Heat         3.64 0.346 15.5     2.91     4.38
# Heat+Water   3.07 0.340 14.6     2.34     3.80
# Water        2.64 0.380 14.6     1.83     3.45
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE   df t.ratio p.value
# Control - Heat           -0.220 0.480 16.6 -0.457  0.9673 
# Control - (Heat+Water)    0.350 0.476 16.1  0.736  0.8812 
# Control - Water           0.783 0.506 16.0  1.550  0.4331 
# Heat - (Heat+Water)       0.570 0.485 15.1  1.175  0.6512 
# Heat - Water              1.003 0.514 15.0  1.951  0.2494 
# (Heat+Water) - Water      0.433 0.510 14.6  0.849  0.8305 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CF15_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Cornflower Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CF15_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CF15_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control        Heat  Heat+Water       Water 
# 0.003344048 0.003684780 0.003071439 0.002638331 
prse
# Control         Heat   Heat+Water        Water 
# 2.736394e-04 3.201829e-04 2.660417e-04 9.523753e-05 
prbarplot=barplot(prmean,   ylim=c(0,0.004), xlab="Treatment",ylab="Mean Cornflower Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 7. 2015 Corn Marigold Seed Number ####

# quick look at the distribution:
plot(density(Seeds$CM15_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glmm1=glmer(CM15_Seeds~Treatment +(1|Plot) + (1|CM15_Day), data=Seeds, family=poisson, na.action = "na.omit")
plot(glmm1)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.13262, p-value = 0.01262
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.6799, p-value < 2.2e-16

# poisson is overdispersed, will try negative binomial:
glmm2=glmer.nb(CM15_Seeds~Treatment +(1|Plot) + (1|CM15_Day), data=Seeds, na.action = "na.omit")
isSingular(glmm2)
# [1] FALSE  # but the variance of Day is small (0.0003048)
plot(glmm2)
fitnegbin <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnegbin)
testUniformity(fitnegbin)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.049204, p-value = 0.8767
testDispersion(fitnegbin)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.9869, p-value = 0.824

# negative binomial is the best fit

summary(glmm2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: Negative Binomial(47.4973)  ( log )
# Formula: CM15_Seeds ~ Treatment + (1 | Plot) + (1 | CM15_Day)
# Data: Seeds
# 
# AIC      BIC   logLik deviance df.resid 
# 1504.4   1525.2   -745.2   1490.4      137 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.40217 -0.60596 -0.06189  0.68928  2.37587 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# Plot     (Intercept) 0.0017165 0.04143 
# CM15_Day (Intercept) 0.0003048 0.01746 
# Number of obs: 144, groups:  Plot, 24; CM15_Day, 2
# 
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)          5.68330    0.03345 169.917  < 2e-16 ***
#   TreatmentHeat       -0.17877    0.04416  -4.048 5.16e-05 ***
#   TreatmentHeat+Water -0.23322    0.04426  -5.269 1.37e-07 ***
#   TreatmentWater      -0.02067    0.04397  -0.470    0.638    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.654              
# TrtmntHt+Wt -0.652  0.494       
# TreatmntWtr -0.657  0.497  0.496

drop1(glmm2, test="Chisq")
# Single term deletions 
# Model:
#   CM15_Seeds ~ Treatment + (1 | Plot) + (1 | CM15_Day)
#           Df    AIC    LRT   Pr(Chi)    
# <none>       1504.4                     
# Treatment  3 1522.1 23.687 2.904e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glmm2, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE  df asymp.LCL asymp.UCL
# Control      5.68 0.0334 Inf      5.60      5.77
# Heat         5.50 0.0337 Inf      5.42      5.59
# Heat+Water   5.45 0.0339 Inf      5.37      5.53
# Water        5.66 0.0335 Inf      5.58      5.75
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat         0.1788 0.0442 Inf  4.048  0.0003 
# Control - Heat+Water   0.2332 0.0443 Inf  5.269  <.0001 
# Control - Water        0.0207 0.0440 Inf  0.470  0.9656 
# Heat - Heat+Water      0.0545 0.0445 Inf  1.224  0.6114 
# Heat - Water          -0.1581 0.0442 Inf -3.579  0.0020 
# Heat+Water - Water    -0.2125 0.0443 Inf -4.799  <.0001 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CM15_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Mean Corn Marigold Seed Number",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CM15_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CM15_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 294.3333   245.9722   233.5556   288.1111
prse
# Control       Heat Heat+Water      Water 
# 8.258810   6.607894   7.579572   6.253733
prbarplot=barplot(prmean,  ylim=c(0,350),xlab="Treatment",ylab="Mean Corn Marigold Seed Number",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 8. 2015 Corn Marigold Average Seed Weight ####


# quick look at the distribution:
plot(density(Seeds$CM15_AWeight, na.rm=TRUE))

mm1 = lmer(CM15_AWeight ~ Treatment +(1|Plot) + (1|CM15_Day), data=Seeds, REML=F, na.action = "na.omit")
isSingular(mm1)
# [1] TRUE 
shapiro.test(residuals(mm1))
# data:  residuals(mm1)
# W = 0.9778, p-value = 0.01931
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)
plot(mm1)

# REML version:
mm1a = lmer(CM15_AWeight ~ Treatment +(1|Plot) + (1|CM15_Day), data=Seeds[-which(is.na(Seeds$CM15_AWeight)),], REML=T, na.action = "na.omit")
isSingular(mm1a)
# [1] TRUE - variance of Day is 0
shapiro.test(residuals(mm1a))
# W = 0.97786, p-value = 0.01963
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

fitg <- simulateResiduals(fittedModel = mm1a)
plot(fitg, asFactor = TRUE)

mm1b = lmer(CM15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CM15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1b))
# W = 0.97786, p-value = 0.01963
mcp.fnc(mm1b)
fitg <- simulateResiduals(fittedModel = mm1b)
plot(fitg)

# this model is OK, worth trying gamma


glmm1=glmer(CM15_AWeight ~ Treatment +(1|Plot) + (1|CM15_Day), data=Seeds[-which(is.na(Seeds$CM15_AWeight)),], family=Gamma, na.action = "na.omit")
isSingular(glmm1)
# TRUE - variance of Day is 0
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.63624, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.13662, p-value < 2.2e-16

# gaussian is the best


summary(mm4b)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CM15_AWeight ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$CM15_AWeight)), ]
# 
# REML criterion at convergence: 148.3
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.96735 -0.71504 -0.08427  0.68432  2.82034 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.008025 0.08958 
# Residual             0.146345 0.38255 
# Number of obs: 144, groups:  Plot, 24
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)           1.6212     0.0735  22.057
# TreatmentHeat        -0.3677     0.1040  -3.537
# TreatmentHeat+Water  -0.3406     0.1040  -3.277
# TreatmentWater       -0.0122     0.1040  -0.117
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.707              
# TrtmntHt+Wt -0.707  0.500       
# TreatmntWtr -0.707  0.500  0.500

drop1(mm1b, test="Chisq")
# Single term deletions
# 
# Model:
#   CM15_AWeight ~ Treatment + (1 | Plot)
#           npar    AIC    LRT   Pr(Chi)    
# <none>         146.37                     
# Treatment    3 158.46 18.093 0.0004208 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(mm1b, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE df lower.CL upper.CL
# Control      1.62 0.0735 20     1.47     1.77
# Heat         1.25 0.0735 20     1.10     1.41
# Heat+Water   1.28 0.0735 20     1.13     1.43
# Water        1.61 0.0735 20     1.46     1.76
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE df t.ratio p.value
# Control - Heat           0.3677 0.104 20  3.537  0.0103 
# Control - (Heat+Water)   0.3406 0.104 20  3.277  0.0182 
# Control - Water          0.0122 0.104 20  0.117  0.9994 
# Heat - (Heat+Water)     -0.0270 0.104 20 -0.260  0.9936 
# Heat - Water            -0.3555 0.104 20 -3.420  0.0133 
# (Heat+Water) - Water    -0.3284 0.104 20 -3.160  0.0234 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Seeds$CM15_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CM15_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CM15_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 1.621209   1.253533   1.280569   1.609014 
prse
# Control       Heat Heat+Water      Water 
# 0.07439168 0.05587479 0.05956930 0.06942786 
prbarplot=barplot(prmean, ylim=c(0,0.002), xlab="Treatment",ylab="Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

# quick graphs of transformed data:
boxplot((Seeds$CM15_AWeight)^(1/3)~Seeds$Treatment,xlab="Treatment",ylab="Cube root Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply((Seeds$CM15_AWeight)^(1/3), Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply((Seeds$CM15_AWeight)^(1/3), Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control        Heat  Heat+Water       Water 
# 0.1165072  0.1069599  0.1077002  0.1163561 
prse
# Control         Heat   Heat+Water        Water 
# 0.001801677 0.001635293 0.001658419 0.001655877
prbarplot=barplot(prmean, ylim=c(0,0.14), xlab="Treatment",ylab="Cube root Mean Corn Marigold Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 9. 2015 Red Dead Nettle Average Seed Weight ####

# quick look at the distribution:
plot(density(Seeds$RDN15_AWeight, na.rm=TRUE))

mm1 = lmer(RDN15_AWeight ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# W = 0.95953, p-value = 0.002063
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(RDN15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$RDN15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# W = 0.9603, p-value = 0.002361
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

fitg <- simulateResiduals(fittedModel = mm1a)
plot(fitg)
testUniformity(fitg)
# D = 0.092727, p-value = 0.3006

# this is OK, worth trying gamma

glmm1=glmer(RDN15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$RDN15_AWeight)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.59983, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.12909, p-value < 2.2e-16


# the gaussian model is the best fit

summary(mm1a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: RDN15_AWeight ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$RDN15_AWeight)), ]
# 
# REML criterion at convergence: -28.8
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.5002 -0.5432  0.1533  0.7388  2.1737 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.005866 0.07659 
# Residual             0.035512 0.18845 
# Number of obs: 110, groups:  Plot, 23
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)          0.81167    0.04649  17.459
# TreatmentHeat       -0.12466    0.06687  -1.864
# TreatmentHeat+Water -0.12306    0.07001  -1.758
# TreatmentWater       0.13061    0.06575   1.987
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.695              
# TrtmntHt+Wt -0.664  0.462       
# TreatmntWtr -0.707  0.492  0.470

drop1(mm1a, test="Chisq")
# Single term deletions
# 
# Model:
#   RDN15_AWeight ~ Treatment + (1 | Plot)
#           npar     AIC    LRT  Pr(Chi)   
# <none>         -34.036                   
# Treatment    3 -24.074 15.962 0.001154 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(mm1a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE   df lower.CL upper.CL
# Control     0.812 0.0465 18.0    0.714    0.909
# Heat        0.687 0.0481 20.1    0.587    0.787
# Heat+Water  0.689 0.0524 19.6    0.579    0.798
# Water       0.942 0.0465 18.0    0.845    1.040
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate     SE   df t.ratio p.value
# Control - Heat           0.1247 0.0669 19.0  1.863  0.2763 
# Control - (Heat+Water)   0.1231 0.0701 18.9  1.756  0.3241 
# Control - Water         -0.1306 0.0657 18.0 -1.987  0.2295 
# Heat - (Heat+Water)     -0.0016 0.0712 19.8 -0.023  1.0000 
# Heat - Water            -0.2553 0.0669 19.0 -3.815  0.0059 
# (Heat+Water) - Water    -0.2537 0.0701 18.9 -3.621  0.0091 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Seeds$RDN15_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Red Dead Nettle Seed Head",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$RDN15_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$RDN15_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.8116667  0.6830247  0.6873188  0.9422778 
prse
# Control       Heat Heat+Water      Water 
# 0.03403148 0.04067062 0.04819385 0.03314925
prbarplot=barplot(prmean, ylim=c(0,0.001), xlab="Treatment",ylab="Mean Seeds per Red Dead Nettle Seed Head",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



#
########################## 10. 2015 Field Speedwell Average Seed Number ####

plot(density(Seeds$SW15_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glmm1=glmer(SW15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Seeds)),], family=poisson, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.098245, p-value = 0.2009
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.2425, p-value = 0.008

# poisson is overdispersed, try neg bin:
glmm2=glmer.nb(SW15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Seeds)),], na.action = "na.omit")
isSingular(glmm2)
# [1] TRUE - the variance of plot is almost 0 (1.771e-11)
plot(glmm2)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitnegbin <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnegbin)
testUniformity(fitnegbin)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.10163, p-value = 0.171
testDispersion(fitnegbin)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.94998, p-value = 0.528

# negative binomial is best, but need to drop Plot:
glm1=glm.nb(SW15_Seeds~Treatment, data=Seeds[-which(is.na(Seeds$SW15_Seeds)),], na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
fitnegbin <- simulateResiduals(fittedModel = glm1, n = 250)
plot(fitnegbin)
testUniformity(fitnegbin)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.09015, p-value = 0.2882
testDispersion(fitnegbin)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.946, p-value = 0.552

summary(glm1)
# Call:
#   glm.nb(formula = SW15_Seeds ~ Treatment, data = Seeds[-which(is.na(Seeds$SW15_Seeds)), 
#                                                         ], na.action = "na.omit", init.theta = 15.66419503, link = log)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -2.85369  -0.65656   0.09406   0.63622   2.15967  
# 
# Coefficients:
#                       Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)          2.43653    0.07102  34.309  < 2e-16 ***
#   TreatmentHeat        0.36578    0.09663   3.785 0.000153 ***
#   TreatmentHeat+Water  0.29132    0.09670   3.013 0.002589 ** 
#   TreatmentWater      -0.07568    0.10157  -0.745 0.456208    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(15.6642) family taken to be 1)
# 
# Null deviance: 159.83  on 118  degrees of freedom
# Residual deviance: 130.18  on 115  degrees of freedom
# AIC: 729.29
# 
# Number of Fisher Scoring iterations: 1 
# 
# Theta:  15.66 
# Std. Err.:  4.83 
# 
# 2 x log-likelihood:  -719.285 

anova(glm1, test="Chisq")
# Model: Negative Binomial(15.6642), link: log 
# Response: SW15_Seeds 
# Terms added sequentially (first to last) 
# Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                        118     159.83                     
# Treatment  3   29.646       115     130.18 9.8821 1.638e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm1, test="Chisq")
# Single term deletions 
# Model:
#   SW15_Seeds ~ Treatment
#           Df Deviance    AIC    LRT  Pr(>Chi)    
# <none>         130.18 727.29                     
# Treatment  3   159.83 750.93 29.646 1.638e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glm1, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE  df asymp.LCL asymp.UCL
# Control      2.44 0.0710 Inf      2.26      2.61
# Heat         2.80 0.0655 Inf      2.64      2.97
# Heat+Water   2.73 0.0656 Inf      2.56      2.89
# Water        2.36 0.0726 Inf      2.18      2.54
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat        -0.3658 0.0966 Inf -3.785  0.0009 
# Control - Heat+Water  -0.2913 0.0967 Inf -3.013  0.0138 
# Control - Water        0.0757 0.1016 Inf  0.745  0.8788 
# Heat - Heat+Water      0.0745 0.0927 Inf  0.803  0.8531 
# Heat - Water           0.4415 0.0978 Inf  4.514  <.0001 
# Heat+Water - Water     0.3670 0.0979 Inf  3.750  0.0010 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Seeds$SW15_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Mean Speedwell Seeds per Head",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$SW15_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$SW15_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 11.43333   16.48276   15.30000   10.60000 
prse
# Control       Heat Heat+Water      Water 
# 0.9810660  0.8576240  0.8725047  0.8088590 
prbarplot=barplot(prmean, ylim=c(0,20), 
                  xlab="Treatment",ylab="Mean Speedwell Seeds per Head",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 11. 2015 Field Speedwell Average Seed Weight ####

plot(density(Seeds$SW15_AWeight, na.rm=TRUE))

mm1 = lmer(SW15_AWeight ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# W = 0.51444, p-value < 2.2e-16
mcp.fnc(mm1)
# not good

# REML = TRUE:
mm1a = lmer(SW15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# W = 0.51444, p-value < 2.2e-16
mcp.fnc(mm1a)
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

fitg <- simulateResiduals(fittedModel = mm1a)
plot(fitg)
testUniformity(fitg)
# D = 0.24114, p-value = 1.952e-06

# this data has very strong right skew, worth trying Gamma

glmm1=glmer(SW15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_AWeight)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.35476, p-value = 1.962e-13
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.71983, p-value = 0.144

# a transformation is probably the best option
# log:
mm2a = lmer(log(SW15_AWeight) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm2a))
# W = 0.95536, p-value = 0.0005809
mcp.fnc(mm2a)
qqnorm(residuals(mm2a))
qqline(residuals(mm2a))
plot(mm2a)
plot(mm2a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm2a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.098353, p-value = 0.1999

# square root:
mm3a = lmer(sqrt(SW15_AWeight) ~ Treatment +(1|Plot) , data=Seeds[-which(is.na(Seeds$SW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm3a))
# W = 0.80651, p-value = 3.231e-11
qqnorm(residuals(mm3a))
qqline(residuals(mm3a))
mcp.fnc(mm3a)
plot(mm3a)
plot(mm3a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm3a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.16793, p-value = 0.002433

# cube root:
mm4a = lmer(SW15_AWeight^(1/3) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm4a))
# W = 0.88271, p-value = 3.123e-08
qqnorm(residuals(mm4a))
qqline(residuals(mm4a))
mcp.fnc(mm4a)
plot(mm4a)
plot(mm4a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm4a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.13234, p-value = 0.03096

# log gaussian looks the best fit

summary(mm2a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: log(SW15_AWeight) ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$SW15_AWeight)), ]
# 
# REML criterion at convergence: 236.1
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.7252 -0.4892  0.1017  0.4713  3.2243 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.05159  0.2271  
# Residual             0.36997  0.6083  
# Number of obs: 119, groups:  Plot, 24
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)         -1.17896    0.14468  -8.149
# TreatmentHeat        0.67562    0.20572   3.284
# TreatmentHeat+Water  0.55855    0.20460   2.730
# TreatmentWater      -0.08162    0.20460  -0.399
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.703              
# TrtmntHt+Wt -0.707  0.497       
# TreatmntWtr -0.707  0.497  0.500

drop1(mm2a, test="Chisq")
# Single term deletions 
# Model:
#   log(SW15_AWeight) ~ Treatment + (1 | Plot)
#           npar    AIC    LRT   Pr(Chi)    
# <none>         239.61                     
# Treatment    3 250.93 17.323 0.0006065 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(mm2a, pairwise~Treatment, adjust="tukey")
# Treatment  lsmean    SE   df lower.CL upper.CL
# Control    -1.179 0.145 19.8   -1.481   -0.877
# Heat       -0.503 0.146 20.6   -0.808   -0.199
# Heat+Water -0.620 0.145 19.8   -0.922   -0.318
# Water      -1.261 0.145 19.8   -1.563   -0.959
# 
# Degrees-of-freedom method: kenward-roger 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE   df t.ratio p.value
# Control - Heat          -0.6756 0.206 20.2 -3.284  0.0178 
# Control - (Heat+Water)  -0.5586 0.205 19.8 -2.730  0.0578 
# Control - Water          0.0816 0.205 19.8  0.399  0.9779 
# Heat - (Heat+Water)      0.1171 0.206 20.2  0.569  0.9402 
# Heat - Water             0.7572 0.206 20.2  3.681  0.0074 
# (Heat+Water) - Water     0.6402 0.205 19.8  3.129  0.0252 
# 
# Degrees-of-freedom method: kenward-roger 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Seeds$SW15_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Speedwell Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$SW15_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$SW15_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control         Heat   Heat+Water        Water 
# 0.0003831201 0.0007441856 0.0006108225 0.0003869231
prse
# Control         Heat   Heat+Water        Water 
# 4.983888e-05 1.544164e-04 6.016080e-05 7.607195e-05 
prbarplot=barplot(prmean, ylim=c(0,0.001), xlab="Treatment",ylab="Mean Speedwell Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

# quick graphs of log data:
boxplot(log(SW15_AWeight) ~ Treatment, data=Seeds, xlab="Treatment",ylab="Mean Speedwell Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(log(Seeds$SW15_AWeight), Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(log(Seeds$SW15_AWeight), Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# -1.1789579 -0.5106635 -0.6204071 -1.2605817
prse
# Control       Heat Heat+Water      Water 
# 0.12422959 0.10401615 0.09513869 0.14200730 
prbarplot=barplot(prmean, ylim=c(-10,0), xlab="Treatment",ylab="Mean Speedwell Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 12. 2015 Common Chickweed Average Seed Number ####


plot(density(Seeds$CW15_Seeds, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glmm1=glmer(CW15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_Seeds)),], family=poisson, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitpois <- simulateResiduals(fittedModel = glmm1, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.12123, p-value = 0.06417
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.82868, p-value = 0.008

# poisson is underdispersed, try negative binomial:
glmm2=glmer.nb(CW15_Seeds~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_Seeds)),], na.action = "na.omit")
# Warning message:
#   In theta.ml(Y, mu, weights = object@resp$weights, limit = limit,  :
#                 iteration limit reached
plot(glmm2)
fitnb <- simulateResiduals(fittedModel = glmm2, n = 250)
plot(fitnb)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.097151, p-value = 0.2194
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.83812, p-value = 0.008

# neg bin is also underdispersed, and there was a warning message.

# try a generalised poisson:
glmm3 = glmmTMB(CW15_Seeds~Treatment +(1|Plot), family=genpois, data=Seeds[-which(is.na(Seeds$CW15_Seeds)),], na.action = "na.omit")
plot(residuals(glmm3, type ="response"))
plot(residuals(glmm3, type ="pearson"))
fitpois <- simulateResiduals(fittedModel = glmm3, n = 250)
plot(fitpois)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.082495, p-value = 0.4034
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.97458, p-value = 0.688

# the generalised poisson looks like the best fit

summary(glmm3)
# Family: genpois  ( log )
# Formula:          CW15_Seeds ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$CW15_Seeds)), ]
# 
# AIC      BIC   logLik deviance df.resid 
# 611.2    627.8   -299.6    599.2      111 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# Plot   (Intercept) 0.007235 0.08506 
# Number of obs: 117, groups:  Plot, 24
# 
# Overdispersion parameter for genpois family (): 0.646 
# 
# Conditional model:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          2.62028    0.05266   49.76  < 2e-16 ***
# TreatmentHeat       -0.27742    0.07919   -3.50  0.00046 ***
# TreatmentHeat+Water -0.16207    0.07624   -2.13  0.03352 *  
# TreatmentWater       0.12643    0.07323    1.73  0.08429 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glmm3, test="Chisq")
# Single term deletions 
# Model:
#   CW15_Seeds ~ Treatment + (1 | Plot)
#           Df    AIC    LRT  Pr(>Chi)    
# <none>       611.18                     
# Treatment  3 625.22 20.035 0.0001669 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glmm3, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE  df lower.CL upper.CL
# Control      2.62 0.0527 111     2.49     2.75
# Heat         2.34 0.0592 111     2.19     2.49
# Heat+Water   2.46 0.0552 111     2.32     2.60
# Water        2.75 0.0509 111     2.62     2.88
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df t.ratio p.value
# Control - Heat          0.277 0.0792 111  3.503  0.0037 
# Control - Heat+Water    0.162 0.0762 111  2.126  0.1513 
# Control - Water        -0.126 0.0732 111 -1.726  0.3150 
# Heat - Heat+Water      -0.115 0.0806 111 -1.430  0.4833 
# Heat - Water           -0.404 0.0781 111 -5.170  <.0001 
# Heat+Water - Water     -0.288 0.0751 111 -3.842  0.0012 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CW15_Seeds~Seeds$Treatment,xlab="Treatment",ylab="Mean Chickweed Seeds per Head",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CW15_Seeds,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CW15_Seeds,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 13.83333   10.29630   11.66667   15.73333 
prse
# Control       Heat Heat+Water      Water 
# 0.4396533  0.7275898  0.6666667  0.3649385 
prbarplot=barplot(prmean, ylim=c(0,20), xlab="Treatment",ylab="Mean Chickweed Seeds per Head",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 13. 2015 Common Chickweed Average Seed Weight ####

plot(density(Seeds$CW15_AWeight, na.rm=TRUE))

mm1 = lmer(CW15_AWeight ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# W = 0.94366, p-value = 9.421e-05
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(CW15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# W = 0.9443, p-value = 0.0001043
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm1a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.10738, p-value = 0.1346


# Gamma:
glmm1=glmer(CW15_AWeight ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_AWeight)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.55598, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.15212, p-value < 2.2e-16

# this data is a bit left-skewed, so an appropriate transformation might be the best option
# square:
mm2a = lmer(CW15_AWeight^2 ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm2a))
# W = 0.9665, p-value = 0.005068
qqnorm(residuals(mm2a))
qqline(residuals(mm2a))
mcp.fnc(mm2a)
plot(mm2a)
plot(mm2a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm2a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.080889, p-value = 0.4282

# cube:
mm3a = lmer(CW15_AWeight^3 ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CW15_AWeight)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm3a))
# W = 0.91131, p-value = 1.02e-06
qqnorm(residuals(mm3a))
qqline(residuals(mm3a))
mcp.fnc(mm3a)
plot(mm3a)
plot(mm3a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm3a, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# D = 0.10414, p-value = 0.158

# squaring the data is the best fit

summary(mm2a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CW15_AWeight^2 ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$CW15_AWeight)), ]
# 
# REML criterion at convergence: -229.3
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.2543 -0.4782  0.0268  0.4835  3.1870 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.001218 0.03490 
# Residual             0.006051 0.07779 
# Number of obs: 117, groups:  Plot, 24
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)         0.155959   0.020117   7.753
# TreatmentHeat       0.010488   0.028864   0.363
# TreatmentHeat+Water 0.029590   0.028450   1.040
# TreatmentWater      0.005411   0.028450   0.190
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.697              
# TrtmntHt+Wt -0.707  0.493       
# TreatmntWtr -0.707  0.493  0.500

drop1(mm2a, test="Chisq")
# Single term deletions
# 
# Model:
#   CW15_AWeight^2 ~ Treatment + (1 | Plot)
#           npar     AIC    LRT Pr(Chi)
# <none>         -241.53               
# Treatment    3 -246.10 1.4323   0.698

lsmeans(mm2a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean     SE   df lower.CL upper.CL
# Control     0.156 0.0201 19.5    0.114    0.198
# Heat        0.166 0.0207 21.7    0.123    0.209
# Heat+Water  0.186 0.0201 19.5    0.144    0.228
# Water       0.161 0.0201 19.5    0.119    0.203
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate     SE   df t.ratio p.value
# Control - Heat         -0.01049 0.0289 20.6 -0.363  0.9831 
# Control - (Heat+Water) -0.02959 0.0284 19.5 -1.040  0.7285 
# Control - Water        -0.00541 0.0284 19.5 -0.190  0.9975 
# Heat - (Heat+Water)    -0.01910 0.0289 20.6 -0.662  0.9102 
# Heat - Water            0.00508 0.0289 20.6  0.176  0.9980 
# (Heat+Water) - Water    0.02418 0.0284 19.5  0.850  0.8300 
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CW15_AWeight~Seeds$Treatment,xlab="Treatment",ylab="Mean Chickweed Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CW15_AWeight,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CW15_AWeight,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.3801819  0.3975702  0.4073777  0.3895470 
prse
# Control       Heat Heat+Water      Water 
# 0.01984511 0.01821313 0.02599244 0.01821666 
prbarplot=barplot(prmean, ylim=c(0,0.0005), xlab="Treatment",ylab="Mean Chickweed Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


# quick graphs of squared data:
boxplot(Seeds$CW15_AWeight2~Seeds$Treatment,xlab="Treatment",ylab="Mean Chickweed Seed Weight",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CW15_AWeight^2 ,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CW15_AWeight^2 ,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.1559593  0.1666868  0.1855492  0.1613704 
prse
# Control       Heat Heat+Water      Water 
# 0.01253659 0.01236205 0.02138152 0.01360724 
prbarplot=barplot(prmean, ylim=c(0,0.00000025), xlab="Treatment",ylab="Mean Chickweed Seed Weight",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)




########################## 14. 2015 Cornflower Nectar Volume ####

plot(density(Seeds$CF15_Nectar, na.rm=TRUE))

mm1 = lmer(CF15_Nectar ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# W = 0.90354, p-value = 0.000157
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(CF15_Nectar ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# W = 0.91399, p-value = 0.0003951
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm1a, n = 250)
plot(fitg, asFactor = TRUE)

# doesn't look good

glmm1=glmer(CF15_Nectar ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.41435, p-value = 5.39e-10
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.42441, p-value < 2.2e-16

# gamma isn't right either, probably need to transform it:
mm2a = lmer(log(CF15_Nectar) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm2a))
# W = 0.97072, p-value = 0.1512
qqnorm(residuals(mm2a))
qqline(residuals(mm2a))
mcp.fnc(mm2a)
plot(mm2a)
plot(mm2a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

# log is pretty good, but maybe a bit strong, try sqrt:
mm3a = lmer(sqrt(CF15_Nectar) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm3a))
# W = 0.98332, p-value = 0.5725
qqnorm(residuals(mm3a))
qqline(residuals(mm3a))
mcp.fnc(mm3a)
plot(mm3a)
plot(mm3a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

# sqrt looks better, worth trying cube root:
mm4a = lmer(CF15_Nectar^(1/3) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm4a))
# W = 0.98931, p-value = 0.8734
qqnorm(residuals(mm4a))
qqline(residuals(mm4a))
mcp.fnc(mm4a)
plot(mm4a)
plot(mm4a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitg <- simulateResiduals(fittedModel = mm4a, n = 250)
plot(fitg, asFactor = TRUE)

glmm2=glmer(CF15_Nectar^(1/3) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$CF15_Nectar)),], family=gaussian(link="log"), na.action = "na.omit")
sapply(c("log", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glm1, family=gaussian(link = x)))
  m_AIC = AIC(update(glm1,  family=gaussian(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#            log   inverse  identity
# [1,]  24.18154  24.18154  24.18154
# [2,] -38.36308 -38.36308 -38.36308
# link function makes no real difference here and is not necessary from an intensive/extensive perspective
# as the data has already been transformed.

# cube root gives the best fit.

summary(mm4a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: CF15_Nectar^(1/3) ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$CF15_Nectar)), ]
# 
# REML criterion at convergence: -80.5
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -2.45673 -0.59658  0.01435  0.50297  2.24764 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.006063 0.07787 
# Residual             0.008878 0.09422 
# Number of obs: 61, groups:  Plot, 18
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)          0.57107    0.04522  12.628
# TreatmentHeat       -0.02854    0.06127  -0.466
# TreatmentHeat+Water -0.05745    0.06542  -0.878
# TreatmentWater       0.06619    0.06393   1.035
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.738              
# TrtmntHt+Wt -0.691  0.510       
# TreatmntWtr -0.707  0.522  0.489

drop1(mm4a, test="Chisq")
# Single term deletions 
# Model:
#   CF15_Nectar^(1/3) ~ Treatment + (1 | Plot)
#           Df     AIC    LRT Pr(Chi)
# <none>       -86.547               
# Treatment  3 -88.008 4.5387  0.2089

lsmeans(mm4a, pairwise~Treatment, adjust="tukey")
# Loading required namespace: pbkrtest
# $lsmeans
# Treatment     lsmean         SE    df  lower.CL  upper.CL
# Control    0.5710707 0.04541974 17.12 0.4752965 0.6668450
# Heat       0.5425316 0.04141205 12.63 0.4527959 0.6322674
# Heat+Water 0.5136224 0.04734133 13.83 0.4119667 0.6152782
# Water      0.6372592 0.04519130 11.71 0.5385258 0.7359926
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate         SE    df t.ratio p.value
# Control - Heat        0.02853912 0.06146471 14.85   0.464  0.9657
# Control - Heat+Water  0.05744830 0.06560605 15.29   0.876  0.8173
# Control - Water      -0.06618844 0.06407188 14.08  -1.033  0.7335
# Heat - Heat+Water     0.02890919 0.06289801 13.29   0.460  0.9665
# Heat - Water         -0.09472755 0.06129610 12.12  -1.545  0.4426
# Heat+Water - Water   -0.12363674 0.06544811 12.76  -1.889  0.2801
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Seeds$CF15_Nectar~Seeds$Treatment,xlab="Treatment",ylab="Cornflower nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CF15_Nectar,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CF15_Nectar,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.2095016  0.1806854  0.1616343  0.2834891
prse
# Control       Heat Heat+Water      Water 
# 0.03072161 0.02935022 0.03336471 0.03728331 
prbarplot=barplot(prmean, ylim=c(0,0.35), xlab="Treatment",ylab="Mean Cornflower nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

# quick graphs of cuberoot data:
boxplot((Seeds$CF15_Nectar)^(1/3)~Seeds$Treatment,xlab="Treatment",ylab="Cornflower nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply((Seeds$CF15_Nectar)^(1/3),Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply((Seeds$CF15_Nectar)^(1/3),Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.5816024  0.5364613  0.5100572  0.6407793 
prse
# Control       Heat Heat+Water      Water 
# 0.02513099 0.02974018 0.03987117 0.02512442 
prbarplot=barplot(prmean, ylim=c(0,0.8), xlab="Treatment",ylab="Mean Cornflower nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)



########################## 15. 2015 Red Dead Nettle Nectar Volume ####

plot(density(Seeds$RDN15_Nectar, na.rm=TRUE))

mm1 = lmer(RDN15_Nectar ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
isSingular(mm1)
# [1] TRUE - variance of Plot is 0
shapiro.test(residuals(mm1))
# 0.92806, p-value = 0.02747
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(RDN15_Nectar ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], REML=T, na.action = "na.omit")
isSingular(mm1a)
# [1] TRUE - variance of Plot is 0
shapiro.test(residuals(mm1a))
# 0.92806, p-value = 0.02747
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)

# Gaussian without Plot:
lm1 = lm(RDN15_Nectar ~ Treatment, data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], na.action = "na.omit")
shapiro.test(residuals(lm1))
# 0.92806, p-value = 0.02747
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
lm1_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Treatment, 
                      Plot=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Plot, 
                      resid=residuals(lm1),
                      fitted=fitted.values(lm1))
ggplot(lm1_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)

# gaussian doesn't look great

# gamma:
glmm1=glmer(RDN15_Nectar~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], family=Gamma, na.action = "na.omit")
isSingular(glmm1)
# [1] TRUE - variance of Plot is 0
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.49882, p-value = 2.406e-08
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.26471, p-value < 2.2e-16

# Gamma without Plot:
glm1=glm(RDN15_Nectar~Treatment, data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], family=Gamma, na.action = "na.omit")
sapply(c("log", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glm1, family=Gamma(link = x)))
  m_AIC = AIC(update(glm1,  family=Gamma(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#            log   inverse  identity
# [1,]  21.68836  21.68836  21.68836
# [2,] -33.37672 -33.37672 -33.37672
# link function made no difference here, so will stick with canonical (inverse)
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Treatment, 
                      Plot=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Plot, 
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.14093, p-value = 0.4671
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.79433, p-value = 0.256

# inverse Gaussian:
glmm2=glmer(RDN15_Nectar~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], family=inverse.gaussian, na.action = "na.omit")
isSingular(glmm2)
# [1] TRUE - variance of Plot is 0
# without Plot:
glm2 = glm(RDN15_Nectar~Treatment, data=Seeds[-which(is.na(Seeds$RDN15_Nectar)),], family=inverse.gaussian, na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Treatment, 
                      Plot=Seeds[-which(is.na(Seeds$RDN15_Nectar)),]$Plot, 
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.17875, p-value = 0.2016
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.57011, p-value = 0.024

# inverse Gaussian is underdispersed, but Gamma looks good

summary(glm1)
# Call:
#   glm(formula = RDN15_Nectar ~ Treatment, family = Gamma, data = Seeds[-which(is.na(Seeds$RDN15_Nectar)), 
#                                                                        ], na.action = "na.omit") 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.17190  -0.40273   0.01423   0.31627   0.92106  
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           2.6610     0.3660   7.271 4.28e-08 ***
# TreatmentHeat         7.0281     2.0967   3.352  0.00218 ** 
# TreatmentHeat+Water   1.2181     0.9940   1.226  0.22992    
# TreatmentWater        0.2194     0.5281   0.415  0.68079    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.2270086)
# 
# Null deviance: 13.2661  on 33  degrees of freedom
# Residual deviance:  8.1632  on 30  degrees of freedom
# AIC: -33.377
# 
# Number of Fisher Scoring iterations: 6

anova(glm1, test="F")
# Analysis of Deviance Table
# Model: Gamma, link: inverse
# Response: RDN15_Nectar
# Terms added sequentially (first to last)
#           Df Deviance Resid. Df Resid. Dev     F    Pr(>F)    
# NULL                         33    13.2661                    
# Treatment  3    5.103        30     8.1632 7.493 0.0006937 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm1, test="F")
# Single term deletions
# Model:
#   RDN15_Nectar ~ Treatment
#           Df Deviance     AIC F value   Pr(>F)   
# <none>         8.1632 -33.377                    
# Treatment  3  13.2661 -16.898  6.2512 0.001998 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glm1, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment    lsmean        SE df asymp.LCL asymp.UCL
# Control    2.661048 0.3660018 NA  1.943698  3.378398
# Heat       9.689104 2.0645228 NA  5.642713 13.735494
# Heat+Water 3.879154 0.9241195 NA  2.067913  5.690395
# Water      2.880414 0.3806312 NA  2.134391  3.626438
# 
# Results are given on the inverse scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate        SE df    z.ratio p.value
# Control - Heat       -7.0280556 2.0967145 NA -3.3519373  0.0045
# Control - Heat+Water -1.2181061 0.9939589 NA -1.2255096  0.6106
# Control - Water      -0.2193662 0.5280506 NA -0.4154264  0.9759
# Heat - Heat+Water     5.8099495 2.2619132 NA  2.5685997  0.0501
# Heat - Water          6.8086894 2.0993176 NA  3.2432869  0.0065
# Heat+Water - Water    0.9987399 0.9994384 NA  0.9993012  0.7498
# 
# P value adjustment: tukey method for comparing a family of 4 estimates


# quick graphs of raw data:
boxplot(Seeds$RDN15_Nectar~Seeds$Treatment,xlab="Treatment",ylab="Red Dead Nettle nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$RDN15_Nectar,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$RDN15_Nectar,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.3757918  0.1032087  0.2577882  0.3471723
prse
# Control       Heat Heat+Water      Water 
# 0.03985445 0.03888758 0.07810586 0.03270765 
prbarplot=barplot(prmean, ylim=c(0,0.5), xlab="Treatment",ylab="Mean Red Dead Nettle nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

(0.1032087/0.3757918)*100
# [1] 27.46433
100 - ((0.1032087/0.3757918)*100)
# [1] 72.53567



########################## 16. 2015 Speedwell Nectar Volume ####

plot(density(Seeds$SW15_Nectar, na.rm=TRUE))

mm1 = lmer(SW15_Nectar ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# 0.97259, p-value = 0.4493
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(SW15_Nectar ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# 0.9639, p-value = 0.24
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
mm1_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Treatment, 
                       Plot=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Plot, 
                       resid=residuals(mm1a),
                       fitted=fitted.values(mm1a))
ggplot(mm1_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitg <- simulateResiduals(fittedModel = mm1a, n = 250)
plot(fitg, asFactor = TRUE)

# gaussian isn't too bad, worth a quick look at others

# gamma:
glmm1=glmer(SW15_Nectar~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Nectar)),], family=Gamma, na.action = "na.omit")
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.52602, p-value = 1.392e-10
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.29498, p-value < 2.2e-16
glmm1_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Treatment, 
                      Plot=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Plot, 
                      resid=residuals(glmm1),
                      fitted=fitted.values(glmm1))
ggplot(glmm1_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)

# transform:
mm2a = lmer(sqrt(SW15_Nectar) ~ Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Nectar)),], REML=T, na.action = "na.omit")
shapiro.test(residuals(mm2a))
# 0.96965, p-value = 0.3663
qqnorm(residuals(mm2a))
qqline(residuals(mm2a))
mcp.fnc(mm2a)
plot(mm2a)
plot(mm2a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
mm2_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Treatment, 
                     Plot=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Plot, 
                     resid=residuals(mm2a),
                     fitted=fitted.values(mm2a))
ggplot(mm2_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
# this only marginally improves the residuals, other transformations made them worse

# gaussian with different (log) link function:
glmm2=glmer(SW15_Nectar~Treatment +(1|Plot), data=Seeds[-which(is.na(Seeds$SW15_Nectar)),], family=gaussian(link="log"), na.action = "na.omit")
shapiro.test(residuals(glmm2))
# 0.97152, p-value = 0.4177
par(mfrow=c(1,1))
qqnorm(residuals(glmm2))
qqline(residuals(glmm2))
mcp.fnc(glmm2)
plot(glmm2)
plot(glmm2, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
m_df <- data.frame(Treatment=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Treatment, 
                     Plot=Seeds[-which(is.na(Seeds$SW15_Nectar)),]$Plot, 
                     resid=residuals(glmm2),
                     fitted=fitted.values(glmm2))
ggplot(m_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitg <- simulateResiduals(fittedModel = glmm2, use.u = T)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.12703, p-value = 0.5143

sapply(c("log", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glmm2, family=gaussian(link = x)))
  m_AIC = AIC(update(glmm2,  family=gaussian(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#             log    inverse   identity
# [1,]   78.02615   77.22424   56.50038
# [2,] -144.05229 -142.44849 -101.00076

# Gaussian with log link on untransformed data gives by far the best model fit.
# a transformative (non-identity) link (such as log) is probably more appropriate here anyway, as the variable 
# of interest (nectar volume) could be considered extensive as I didn't control for flower/nectary size.

summary(glmm2)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: gaussian  ( log )
# Formula: SW15_Nectar ~ Treatment + (1 | Plot)
# Data: Seeds[-which(is.na(Seeds$SW15_Nectar)), ]
# 
# AIC      BIC   logLik deviance df.resid 
# -144.1   -134.1     78.0   -156.1       33 
# 
# Scaled residuals: 
#   Min       1Q   Median       3Q      Max 
# -1.97721 -0.43021 -0.01734  0.48019  2.22525 
# 
# Random effects:
#   Groups   Name        Variance  Std.Dev.
# Plot     (Intercept) 0.0807168 0.28411 
# Residual             0.0009624 0.03102 
# Number of obs: 39, groups:  Plot, 22
# 
# Fixed effects:
#                     Estimate Std. Error t value Pr(>|z|)    
# (Intercept)         -2.13741    0.23584  -9.063  < 2e-16 ***
# TreatmentHeat       -0.90914    0.33162  -2.741  0.00612 ** 
# TreatmentHeat+Water -1.15107    0.36431  -3.160  0.00158 ** 
# TreatmentWater      -0.01666    0.30271  -0.055  0.95612    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.685              
# TrtmntHt+Wt -0.635  0.438       
# TreatmntWtr -0.708  0.498  0.456

drop1(glmm2, test="Chisq")
# Single term deletions
# Model:
#   SW15_Nectar ~ Treatment + (1 | Plot)
#           Df     AIC    LRT Pr(Chi)  
# <none>       -144.05                 
# Treatment  3 -139.56 10.497 0.01478 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(glmm2, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment  lsmean    SE  df asymp.LCL asymp.UCL
# Control     -2.14 0.236 Inf     -2.72     -1.55
# Heat        -3.05 0.242 Inf     -3.65     -2.44
# Heat+Water  -3.29 0.282 Inf     -3.99     -2.59
# Water       -2.15 0.215 Inf     -2.69     -1.62
# 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate    SE  df z.ratio p.value
# Control - Heat         0.9091 0.332 Inf  2.741  0.0311 
# Control - Heat+Water   1.1511 0.364 Inf  3.160  0.0086 
# Control - Water        0.0167 0.303 Inf  0.055  0.9999 
# Heat - Heat+Water      0.2419 0.370 Inf  0.654  0.9142 
# Heat - Water          -0.8925 0.319 Inf -2.799  0.0264 
# Heat+Water - Water    -1.1344 0.352 Inf -3.224  0.0069 
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Seeds$SW15_Nectar~Seeds$Treatment,xlab="Treatment",ylab="Speedwell nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$SW15_Nectar,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$SW15_Nectar,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 0.14256231 0.05035047 0.03792835 0.12689040 
prse
# Control        Heat  Heat+Water       Water 
# 0.016758479 0.009261401 0.005743284 0.016703763
prbarplot=barplot(prmean, ylim=c(0,0.2), xlab="Treatment",ylab="Mean Speedwell nectar volume (microL)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

(0.05035047/0.14256231)*100
# [1] 35.31822
100 - ((0.05035047/0.14256231)*100)
# [1] 64.68178




#
########################## 17. Corn Marigold Flower Diameter ####

plot(density(Seeds$CM15_DiscDiameter, na.rm=TRUE))

mm1 = lmer(CM15_DiscDiameter ~ Treatment +(1|Plot), data=Seeds, REML=F, na.action = "na.omit")
shapiro.test(residuals(mm1))
# 0.9605, p-value = 3.588e-06
qqnorm(residuals(mm1))
qqline(residuals(mm1))
mcp.fnc(mm1)

# REML = TRUE:
mm1a = lmer(CM15_DiscDiameter ~ Treatment +(1|Plot), data=Seeds, REML=T, na.action = "na.omit")
shapiro.test(residuals(mm1a))
# 0.96045, p-value = 3.539e-06
qqnorm(residuals(mm1a))
qqline(residuals(mm1a))
mcp.fnc(mm1a)
plot(mm1a)
plot(mm1a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
mm1_df <- data.frame(Treatment=Seeds$Treatment, 
                     Plot=Seeds$Plot, 
                     resid=residuals(mm1a),
                     fitted=fitted.values(mm1a))
ggplot(mm1_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
ggplot(mm1_df,aes(x=Plot,y=resid))+geom_point()+geom_boxplot(aes(group=Plot),alpha = 0)
# it looks like it would be a perfect fit for Gaussian but for the two outliers
fitg <- simulateResiduals(fittedModel = mm1a)
plot(fitg, asFactor = TRUE)
testUniformity(fitg) 
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.069511, p-value = 0.1965

# worth trying Gamma

# Gamma:
glmm1=glmer(CM15_DiscDiameter~Treatment +(1|Plot), data=Seeds, family=Gamma, na.action = "na.omit")
# Warning message:
#   In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                  Model is nearly unidentifiable: very large eigenvalue
#                - Rescale variables?
plot(glmm1)
plot(glmm1, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
fitgamma <- simulateResiduals(fittedModel = glmm1, use.u = T)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.75034, p-value < 2.2e-16
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.044534, p-value < 2.2e-16

# try transformations:
mm2a = lmer(log(CM15_DiscDiameter) ~ Treatment +(1|Plot), data=Seeds, REML=T, na.action = "na.omit")
shapiro.test(residuals(mm2a))
# 0.97788, p-value = 0.000838
qqnorm(residuals(mm2a))
qqline(residuals(mm2a))
mcp.fnc(mm2a)
plot(mm2a)
plot(mm2a, resid(., scaled=TRUE) ~ fitted(.) | Treatment , abline = 0)
mm2_df <- data.frame(Treatment=Seeds$Treatment, 
                     Plot=Seeds$Plot, 
                     resid=residuals(mm2a),
                     fitted=fitted.values(mm2a))
ggplot(mm2_df,aes(x=Treatment,y=resid))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
ggplot(mm2_df,aes(x=Plot,y=resid))+geom_point()+geom_boxplot(aes(group=Plot),alpha = 0)
fitg <- simulateResiduals(fittedModel = mm2a)
# Model family was recognized or set as continuous, but duplicate values were detected in the response. Consider if you are fitting an appropriate model.
plot(fitg)
testUniformity(fitg) 
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.052163, p-value = 0.531

# log Gaussian is the best

summary(mm2a)
# Linear mixed model fit by REML ['lmerMod']
# Formula: log(CM15_DiscDiameter) ~ Treatment + (1 | Plot)
# Data: Seeds
# 
# REML criterion at convergence: -316
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.0238 -0.5814  0.0704  0.5523  4.5528 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Plot     (Intercept) 0.005444 0.07378 
# Residual             0.012415 0.11142 
# Number of obs: 240, groups:  Plot, 24
# 
# Fixed effects:
#                     Estimate Std. Error t value
# (Intercept)          2.64898    0.03338  79.360
# TreatmentHeat       -0.08988    0.04721  -1.904
# TreatmentHeat+Water -0.07437    0.04721  -1.575
# TreatmentWater       0.06480    0.04721   1.373
# 
# Correlation of Fixed Effects:
#   (Intr) TrtmnH TrtH+W
# TreatmentHt -0.707              
# TrtmntHt+Wt -0.707  0.500       
# TreatmntWtr -0.707  0.500  0.500

drop1(mm2a, test="Chisq")
# Single term deletions 
# Model:
#   log(CM15_DiscDiameter) ~ Treatment + (1 | Plot)
#           Df     AIC    LRT  Pr(Chi)   
# <none>       -324.21                   
# Treatment  3 -317.64 12.568 0.005671 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

lsmeans(mm2a, pairwise~Treatment, adjust="tukey")
# $lsmeans
# Treatment    lsmean         SE df lower.CL upper.CL
# Control    2.648985 0.03337931 20 2.579357 2.718613
# Heat       2.559103 0.03337931 20 2.489475 2.628731
# Heat+Water 2.574616 0.03337931 20 2.504988 2.644244
# Water      2.713782 0.03337931 20 2.644154 2.783410
# 
# Results are given on the log scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate         SE df t.ratio p.value
# Control - Heat        0.08988219 0.04720548 20   1.904  0.2579
# Control - Heat+Water  0.07436885 0.04720548 20   1.575  0.4142
# Control - Water      -0.06479687 0.04720548 20  -1.373  0.5300
# Heat - Heat+Water    -0.01551335 0.04720548 20  -0.329  0.9874
# Heat - Water         -0.15467906 0.04720548 20  -3.277  0.0182
# Heat+Water - Water   -0.13916571 0.04720548 20  -2.948  0.0367
# 
# Results are given on the log scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates


# quick graphs of raw data:
boxplot(Seeds$CM15_DiscDiameter~Seeds$Treatment,xlab="Treatment",ylab="Corn Marigold flower disc diameter (mm)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Seeds$CM15_DiscDiameter,Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(Seeds$CM15_DiscDiameter,Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 14.24500   13.01167   13.29000   15.19333 
prse
# Control       Heat Heat+Water      Water 
# 0.2269974  0.1948358  0.2802108  0.2364625
prbarplot=barplot(prmean, ylim=c(0,20), xlab="Treatment",ylab="Mean Corn Marigold flower disc diameter (mm)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


# quick graphs of log data:
boxplot(log(Seeds$CM15_DiscDiameter)~Seeds$Treatment,xlab="Treatment",ylab="Corn Marigold flower disc diameter (mm)",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(log(Seeds$CM15_DiscDiameter),Seeds$Treatment, na.rm=TRUE, mean)
prse=tapply(log(Seeds$CM15_DiscDiameter),Seeds$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 2.648985   2.559103   2.574616   2.713782 
prse
# Control       Heat Heat+Water      Water 
# 0.01584627 0.01521366 0.02032300 0.01547907
prbarplot=barplot(prmean, ylim=c(0,3), xlab="Treatment",ylab="Mean Corn Marigold flower disc diameter (mm)",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)

