.libPaths(c("E:\\Rlib", .libPaths()))

library(emmeans)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(plotrix)


######################### Functions


######################### Datasets


# 1. Total Floral Abundance
# 2. Floral Richness
# 3. Visitor Abundance (no. of visits) 
# 4. Visitor Richness (extrapolated)
# 5. Visits Per Flower (all species)
# 6. Visits Per C. cyanus Flower
# 7. Visits Per G. segetum Flower
# 8. Visits Per L. purpureum Flower
# 9. Visits Per V. persica Flower
# 10. Visits Per S. media Flower
# 11. Diet Breadth
# 12. Weighted Connectance
# 13. Interaction Evenness
# 14. Generality
# 15. Vulnerability



# load the data file:
Flowers = read.csv("Data/AllData_PlotLevel_BothYears.csv",header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
str(Flowers)
# change plot to factor:
Flowers$Plot = as.factor(Flowers$Plot)
# change Treatment to factor:
Flowers$Treatment = as.factor(Flowers$Treatment)
# Year to a factor:
Flowers$Year = as.factor(Flowers$Year)
# some of the columns have been classed as character, need to set all relevant ones to numeric:
Flowers[, c(23, 29, 32)] = lapply(Flowers[, c(23, 29, 32)], as.numeric)



########################## 1. Total Floral Abundance ####

# quick look at the distribution:
plot(density(Flowers$TotalFlowers, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glm1=glm(TotalFlowers~Treatment*Year, data=Flowers, family=poisson)
glm1$deviance / glm1$df.residual
# [1] 162.3002
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
fitpois <- simulateResiduals(fittedModel = glm1, n = 250)
plot(fitpois, asFactor = TRUE)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.41667, p-value = 1.156e-07
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 11.451, p-value < 2.2e-16

# poisson was very overdispersed, have a look at negative binomial:
glm2=glm.nb(TotalFlowers~Treatment*Year, data=Flowers)
glm2$deviance / glm2$df.residual
# [1] 1.22403
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
fitnb <- simulateResiduals(fittedModel = glm2, n = 250)
plot(fitnb, asFactor = TRUE)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.1263, p-value = 0.3952
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.88878, p-value = 0.4

# negative binomial looks good

summary(glm2)
# Call:
#   glm.nb(formula = TotalFlowers ~ Treatment * Year, data = Flowers, 
#          init.theta = 8.613501995, link = log)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -2.33958  -0.82455   0.05118   0.57065   1.97205  
# 
# Coefficients:
#                           Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 7.3497     0.1395  52.691  < 2e-16 ***
# TreatmentHeat              -0.6304     0.1975  -3.192  0.00141 ** 
# TreatmentHeat+Water        -0.5213     0.1975  -2.640  0.00829 ** 
# TreatmentWater             -0.2213     0.1973  -1.121  0.26210    
# Year2                       0.1126     0.1972   0.571  0.56823    
# TreatmentHeat:Year2         0.2841     0.2792   1.018  0.30887    
# TreatmentHeat+Water:Year2   0.1397     0.2791   0.500  0.61680    
# TreatmentWater:Year2        0.3683     0.2790   1.320  0.18676    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(8.6135) family taken to be 1)
# 
# Null deviance: 81.705  on 47  degrees of freedom
# Residual deviance: 48.961  on 40  degrees of freedom
# AIC: 734.92
# 
# Number of Fisher Scoring iterations: 1 
# 
# Theta:  8.61 
# Std. Err.:  1.74 
# 
# 2 x log-likelihood:  -716.917 

anova(glm2, test="Chisq")
# Analysis of Deviance Table 
# Model: Negative Binomial(8.6135), link: log 
# Response: TotalFlowers 
# Terms added sequentially (first to last) 
#                Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                              47     81.705              
# Treatment       3  20.9014        44     60.804 0.0001104 ***
# Year            1   9.8151        43     50.989 0.0017309 ** 
# Treatment:Year  3   2.0273        40     48.961 0.5667564    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Warning message:
#   In anova.negbin(glm2, test = "Chisq") :
#   tests made without re-estimating 'theta'


# the interaction term is not significant, so I will re-run without it:
glm2b=glm.nb(TotalFlowers~Treatment+Year, data=Flowers, link = "log")
# worth checking the link function doesn't improve model fit before assessing main effects:
poisson_links = c("log", "sqrt", "identity")
logLik(update(glm2b, link = poisson_links[1]))
# 'log Lik.' -359.4521 (df=6)
logLik(update(glm2b, link = poisson_links[2]))
# 'log Lik.' -359.4269 (df=6)
logLik(update(glm2b, link = poisson_links[3]))
# 'log Lik.' -359.5111 (df=6)
AIC(update(glm2b, link = poisson_links[1]))
# [1] 730.9043
AIC(update(glm2b, link = poisson_links[2]))
# [1] 730.8539
AIC(update(glm2b, link = poisson_links[3]))
# [1] 731.0221
# these are all extremely similar, plots also showed very little difference. Going to stick with canonical (log).
par(mfrow=c(2,2))
plot(glm2b)
par(mfrow=c(1,1))
m_df <- data.frame(Treatment=Flowers$Treatment, 
                   Year=Flowers$Year,
                   TrYe=paste0(Flowers$Treatment, Flowers$Year),
                   resid=residuals(glm2b),
                   fitted=fitted.values(glm2b))
ggplot(m_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(m_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitnb <- simulateResiduals(fittedModel = glm2b, n = 250)
plot(fitnb, asFactor = TRUE)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.11682, p-value = 0.4927
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.90088, p-value = 0.528

AIC(glm2, glm2b)
# df      AIC
# glm2   9 734.9173
# glm2b  6 730.9043

anova(glm2, glm2b)
# Likelihood ratio tests of Negative Binomial Models 
# Response: TotalFlowers
#              Model    theta Resid. df    2 x log-lik.   Test    df LR stat.   Pr(Chi)
# 1 Treatment + Year 8.273710        43       -718.9043                                
# 2 Treatment * Year 8.613502        40       -716.9173 1 vs 2     3 1.987008 0.5751077

summary(glm2b)
# Call:
#   glm.nb(formula = TotalFlowers ~ Treatment + Year, data = Flowers, 
#          init.theta = 8.273710028, link = log)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -2.48801  -0.66206   0.07331   0.45759   1.93295  
# 
# Coefficients:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          7.25547    0.11252  64.479  < 2e-16 ***
# TreatmentHeat       -0.49224    0.14240  -3.457 0.000547 ***
# TreatmentHeat+Water -0.45589    0.14239  -3.202 0.001366 ** 
# TreatmentWater      -0.03828    0.14230  -0.269 0.787893    
# Year2                0.31057    0.10069   3.084 0.002040 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(8.2737) family taken to be 1)
# 
# Null deviance: 78.506  on 47  degrees of freedom
# Residual deviance: 48.993  on 43  degrees of freedom
# AIC: 730.9
# 
# Number of Fisher Scoring iterations: 1 
# 
# Theta:  8.27 
# Std. Err.:  1.67 
# 
# 2 x log-likelihood:  -718.904

# Type II LRT:
drop1(glm2b, test="Chisq")
# Single term deletions
# 
# Model:
#   TotalFlowers ~ Treatment + Year
#           Df Deviance    AIC     LRT  Pr(>Chi)    
# <none>         48.993 728.90                      
# Treatment  3   69.371 743.28 20.3775 0.0001418 ***
# Year       1   58.424 736.33  9.4305 0.0021340 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Type I LRT:
anova(glm2b, test="Chisq")
# Analysis of Deviance Table 
# Model: Negative Binomial(8.2737), link: log 
# Response: TotalFlowers 
# Terms added sequentially (first to last) 
# 
#           Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                         47     78.506              
# Treatment  3  20.0820        44     58.424 0.0001632 ***
# Year       1   9.4305        43     48.993 0.0021340 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm2b, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean    SE  df asymp.LCL asymp.UCL
# Control      7.41 0.101 Inf      7.16      7.66
# Heat         6.92 0.101 Inf      6.67      7.17
# Heat+Water   6.95 0.101 Inf      6.70      7.21
# Water        7.37 0.101 Inf      7.12      7.62
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate    SE  df z.ratio p.value
# Control - Heat         0.4922 0.142 Inf  3.457  0.0031 
# Control - Heat+Water   0.4559 0.142 Inf  3.202  0.0075 
# Control - Water        0.0383 0.142 Inf  0.269  0.9932 
# Heat - Heat+Water     -0.0363 0.143 Inf -0.255  0.9942 
# Heat - Water          -0.4540 0.142 Inf -3.188  0.0078 
# Heat+Water - Water    -0.4176 0.142 Inf -2.933  0.0177 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$TotalFlowers~Flowers$Treatment,xlab="Treatment",ylab="Floral abundance",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$TotalFlowers,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$TotalFlowers,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 1648.333   1029.750   1056.167   1631.750 
prse
#   Control       Heat Heat+Water      Water 
# 143.2510   101.1571   165.4637   151.5542 
prbarplot=barplot(prmean, ylim=c(0,2000), xlab="Treatment",ylab="Mean floral abundance",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(TotalFlowers~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Floral abundance", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$TotalFlowers, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$TotalFlowers,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 1555.667  828.1667   923.6667 1246.833
# 2 1741.000 1231.3333  1188.6667 2016.667
prse
#    Control     Heat Heat+Water    Water 
# 1 197.3933 123.4942   190.5532 133.4744
# 2 218.8441 116.2711   277.7251 154.8627
prbarplot=barplot(prmean, beside=T, ylim=c(0,2500), xlab="Treatment",ylab="Mean Floral Abundance", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 2500, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))


# Heat:
100 - ((1029.750/1648.333)*100)
# [1] 37.52779

# Heat+Water
100 - ((1056.167/1648.333)*100)
# [1] 35.92514


#
########################## 2. Floral Richness ####

# quick look at the distribution:
plot(density(Flowers$FRichness, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glm1=glm(FRichness~Treatment*Year, data=Flowers, family=poisson)
glm1$deviance / glm1$df.residual
# [1] 0.2936171
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
fitpois <- simulateResiduals(fittedModel = glm1, n = 250)
plot(fitpois, asFactor = TRUE)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.21947, p-value = 0.01644
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.47589, p-value < 2.2e-16

# negative binomial:
glm2=glm.nb(FRichness~Treatment*Year, data=Flowers)
# Warning messages:
#   1: In theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace >  :
#                    iteration limit reached
glm2$deviance / glm2$df.residual
# [1] 0.2936129
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
fitnb <- simulateResiduals(fittedModel = glm2, n = 250)
plot(fitnb, asFactor = TRUE)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.21991, p-value = 0.01613
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.47536, p-value < 2.2e-16

# poisson and negative binomial are underdispersed, will try generalised poisson:
glm3 = glmmTMB(FRichness~Treatment*Year, data=Flowers, family=genpois, na.action = "na.omit")
plot(residuals(glm3, type ="response"))
plot(residuals(glm3, type ="pearson"))
glm3_df <- data.frame(Treatment=Flowers$Treatment, 
                     Year=Flowers$Year,
                     TrYe=paste0(Flowers$Treatment, Flowers$Year),
                     resid=residuals(glm3, type ="pearson"))
ggplot(glm4_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm4_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgpois <- simulateResiduals(fittedModel = glm3, n = 250)
plot(fitgpois, asFactor = TRUE)
testUniformity(fitgpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.073988, p-value = 0.9378
testDispersion(fitgpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.94272, p-value = 0.584

# generalised poisson is the best fit

summary(glm3)
# Family: genpois  ( log )
# Formula:          FRichness ~ Treatment * Year
# Data: Flowers
# 
# AIC      BIC   logLik deviance df.resid 
# 204.3    221.1    -93.2    186.3       39 
# 
# 
# Overdispersion parameter for genpois family (): 0.26 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                2.71563    0.05339   50.86  < 2e-16 ***
# TreatmentHeat             -0.13870    0.07724   -1.80   0.0725 .  
# TreatmentHeat+Water       -0.10362    0.07697   -1.35   0.1782    
# TreatmentWater            -0.09319    0.07735   -1.20   0.2283    
# Year2                     -0.50768    0.08625   -5.89 3.95e-09 ***
# TreatmentHeat:Year2        0.03290    0.12219    0.27   0.7878    
# TreatmentHeat+Water:Year2  0.17218    0.12099    1.42   0.1547    
# TreatmentWater:Year2       0.16870    0.12214    1.38   0.1672    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm3, test="Chisq")
# Single term deletions
# Model:
#   FRichness ~ Treatment * Year
#                Df    AIC    LRT Pr(>Chi)
# <none>            204.31                
# Treatment:Year  3 201.51 3.2052   0.3611


# the interaction is not significant, so I will drop it from the model
# generalised poisson without interaction term:
glm4 = glmmTMB(FRichness~Treatment+Year, data=Flowers, family=genpois, na.action = "na.omit")
sapply(c("log", "sqrt", "identity"), function(x) {
  m_loglik = logLik(update(glm4, family=poisson(link = x)))
  m_AIC = AIC(update(glm4,  family=poisson(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#            log      sqrt  identity
# [1,] -108.7030 -108.6433 -108.5982
# [2,]  227.4059  227.2865  227.1964
# these are extremely similar, diagnostic plots showed link function makes little difference. Will stick
# with canonical (log).
plot(residuals(glm4, type ="response"))
plot(residuals(glm4, type ="pearson"))
glm4_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm4, type ="pearson"))
ggplot(glm4_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm4_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgpois <- simulateResiduals(fittedModel = glm4, n = 250)
plot(fitgpois, asFactor = TRUE)
testUniformity(fitgpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.068105, p-value = 0.9681
testDispersion(fitgpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 0.96009, p-value = 0.712

summary(glm4)
# Family: genpois  ( log )
# Formula:          FRichness ~ Treatment + Year
# Data: Flowers
# 
# AIC      BIC   logLik deviance df.resid 
# 201.5    212.7    -94.8    189.5       42 
# 
# 
# Overdispersion parameter for genpois family (): 0.274 
# 
# Conditional model:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          2.68033    0.04599   58.29   <2e-16 ***
# TreatmentHeat       -0.13055    0.06226   -2.10    0.036 *  
# TreatmentHeat+Water -0.03543    0.06114   -0.58    0.562    
# TreatmentWater      -0.02603    0.06151   -0.42    0.672    
# Year2               -0.41409    0.04459   -9.29   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm4, test="Chisq")
# Single term deletions 
# Model:
#   FRichness ~ Treatment + Year
#           Df    AIC    LRT  Pr(>Chi)    
# <none>       201.51                     
# Treatment  3 200.45  4.934    0.1767    
# Year       1 250.58 51.069 8.916e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm4, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE df lower.CL upper.CL
# Control      2.47 0.0436 42     2.36     2.59
# Heat         2.34 0.0457 42     2.22     2.46
# Heat+Water   2.44 0.0435 42     2.32     2.55
# Water        2.45 0.0439 42     2.33     2.56
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE df t.ratio p.value
# Control - Heat         0.1305 0.0623 42  2.097  0.1708 
# Control - Heat+Water   0.0354 0.0611 42  0.580  0.9377 
# Control - Water        0.0260 0.0615 42  0.423  0.9742 
# Heat - Heat+Water     -0.0951 0.0619 42 -1.536  0.4256 
# Heat - Water          -0.1045 0.0629 42 -1.661  0.3567 
# Heat+Water - Water    -0.0094 0.0616 42 -0.153  0.9987 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$FRichness~Flowers$Treatment,xlab="Treatment",ylab="Floral abundance",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$FRichness,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$FRichness,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 12.25000   10.33333   11.66667   12.00000 
prse
#   Control       Heat Heat+Water      Water 
# 0.9856931  1.0683699  0.7913676  0.6396021 
prbarplot=barplot(prmean, ylim=c(0,15), xlab="Treatment",ylab="Mean floral abundance",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(FRichness~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Floral abundance", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$FRichness, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$FRichness,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 15.333333 13.000000  13.666667    14
# 2  9.166667  7.666667   9.666667    10
prse
#    Control     Heat Heat+Water    Water 
# 1 0.3333333 1.000000  0.7149204 0.2581989
# 2 0.6009252 1.085255  0.8027730 0.3651484
prbarplot=barplot(prmean, beside=T, ylim=c(0,20), xlab="Treatment",ylab="Mean Floral Abundance", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 2500, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



########################## 3. Visitor Abundance (no. of visits) ####

# quick look at the distribution:
plot(density(Flowers$Visits, na.rm=TRUE))

# this is count data so I will start with poisson

# poisson:
glm1=glm(Visits ~ Treatment*Year, data=Flowers, family=poisson, na.action = "na.omit")
glm1$deviance / glm1$df.residual
# [1] 4.167082 # looks overdispersed
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
fitpois <- simulateResiduals(fittedModel = glm1, n = 250)
plot(fitpois, asFactor = TRUE)
testUniformity(fitpois)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.24477, p-value = 0.006354
testDispersion(fitpois)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.8944, p-value < 2.2e-16

# poisson is overdispersed, try negative binomial:
glm2=glm.nb(Visits ~ Treatment*Year, data=Flowers)
glm2$deviance / glm2$df.residual
# [1] 1.22403
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
fitnb <- simulateResiduals(fittedModel = glm2, n = 250)
plot(fitnb, asFactor = TRUE)
testUniformity(fitnb)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.09877, p-value = 0.7001
testDispersion(fitnb)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated 
# data:  simulationOutput
# ratioObsSim = 1.0074, p-value = 0.888

# negative binomial looks good

summary(glm2)
# glm.nb(formula = Visits ~ Treatment * Year, data = Flowers, init.theta = 31.96239346, 
#        link = log)
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -2.14787  -0.79770  -0.08525   0.74023   2.19268  
# 
# Coefficients:
#                             Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                4.449e+00  8.464e-02  52.559  < 2e-16 ***
# TreatmentHeat             -3.625e-01  1.232e-01  -2.943  0.00325 ** 
# TreatmentHeat+Water       -3.376e-01  1.229e-01  -2.747  0.00602 ** 
# TreatmentWater            -1.288e-01  1.208e-01  -1.066  0.28634    
# Year2                     -2.125e-15  1.197e-01   0.000  1.00000    
# TreatmentHeat:Year2        4.956e-01  1.711e-01   2.897  0.00377 ** 
# TreatmentHeat+Water:Year2  4.199e-01  1.711e-01   2.454  0.01413 *  
# TreatmentWater:Year2       1.782e-01  1.698e-01   1.050  0.29381    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Negative Binomial(31.9624) family taken to be 1)
# 
# Null deviance: 81.072  on 47  degrees of freedom
# Residual deviance: 48.769  on 40  degrees of freedom
# AIC: 423.97
# 
# Number of Fisher Scoring iterations: 1 
# 
# Theta:  31.96 
# Std. Err.:  9.25 
# 
# 2 x log-likelihood:  -405.974 

anova(glm2, test="Chisq")
# Analysis of Deviance Table# 
# Model: Negative Binomial(31.9624), link: log# 
# Response: Visits# 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                              47     81.072              
# Treatment       3   1.8823        44     79.190   0.59718    
# Year            1  19.8443        43     59.345 8.401e-06 ***
# Treatment:Year  3  10.5760        40     48.769   0.01425 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm2, pairwise~Treatment, adjust="tukey")
# NOTE: Results may be misleading due to involvement in interactions
# $emmeans
# Treatment  emmean     SE  df asymp.LCL asymp.UCL
# Control      4.45 0.0598 Inf      4.30      4.60
# Heat         4.33 0.0611 Inf      4.18      4.49
# Heat+Water   4.32 0.0611 Inf      4.17      4.47
# Water        4.41 0.0602 Inf      4.26      4.56
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat         0.1147 0.0855 Inf  1.342  0.5363 
# Control - Heat+Water   0.1277 0.0856 Inf  1.492  0.4422 
# Control - Water        0.0397 0.0849 Inf  0.467  0.9662 
# Heat - Heat+Water      0.0129 0.0864 Inf  0.150  0.9988 
# Heat - Water          -0.0751 0.0858 Inf -0.875  0.8178 
# Heat+Water - Water    -0.0880 0.0858 Inf -1.025  0.7346 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Flowers$Visits~Flowers$Treatment,xlab="Treatment",ylab="Visitor abundance",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$Visits,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$Visits,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 85.50000   78.58333   76.91667   82.50000
prse
#  Control       Heat Heat+Water      Water 
# 4.635992   8.013680   7.426874   4.622901 
prbarplot=barplot(prmean, ylim=c(0,120), xlab="Treatment",ylab="Mean Visitor abundance",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(Visits~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Visitor abundance", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$Visits, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$Visits,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#   Control     Heat Heat+Water    Water 
# 1    85.5 59.50000   61.00000 75.16667
# 2    85.5 97.66667   92.83333 89.83333
prse
#    Control     Heat Heat+Water    Water
# 1 5.155903 8.192476   5.802298 4.847107
# 2 8.245201 8.353309  10.377593 7.001984
prbarplot=barplot(prmean, beside=T, ylim=c(0,120), xlab="Treatment",ylab="Mean Visitor abundance", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 120, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



########################## 4. Visitor Richness (extrapolated) ####

# quick look at the distribution:
plot(density(Flowers$Chao, na.rm=TRUE))

# this data is non-integer, so poisson and neg bin are not appropriate. Will start with gaussian:

lm1=lm(Chao ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.90998, p-value = 0.001342
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)
testUniformity(fitg)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.15483, p-value = 0.2

# try gamma:
glm1=glm(Chao ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
# link function made no difference here
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.082, p-value = 0.9035
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.1907, p-value = 0.408

# quick look at inverse gaussian:
glm2=glm(Chao ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
# link function made no difference here
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.094833, p-value = 0.7811
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.85733, p-value = 0.872

# gamma looks the best

summary(glm1)
# Call:
#   glm(formula = Chao ~ Treatment * Year, family = Gamma(link = "log"), 
#       data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.97785  -0.34215  -0.08422   0.13273   0.96482  
# 
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                3.586447   0.203157  17.654   <2e-16 ***
# TreatmentHeat              0.172120   0.287307   0.599    0.552    
# TreatmentHeat+Water       -0.157211   0.287307  -0.547    0.587    
# TreatmentWater             0.105794   0.287307   0.368    0.715    
# Year2                      0.001268   0.287307   0.004    0.997    
# TreatmentHeat:Year2       -0.325079   0.406313  -0.800    0.428    
# TreatmentHeat+Water:Year2  0.273294   0.406313   0.673    0.505    
# TreatmentWater:Year2       0.099824   0.406313   0.246    0.807    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.2476358)
# 
# Null deviance: 10.0781  on 47  degrees of freedom
# Residual deviance:  9.2855  on 40  degrees of freedom
# AIC: 415.46
# 
# Number of Fisher Scoring iterations: 5

anova(glm1, test="F")
# Analysis of Deviance Table
# Model: Gamma, link: log 
# Response: Chao 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                              47    10.0781              
# Treatment       3  0.22331        44     9.8548 0.3006 0.8248
# Year            1  0.00222        43     9.8526 0.0090 0.9250
# Treatment:Year  3  0.56707        40     9.2855 0.7633 0.5213

# the interaction is not significant, so will rerun without it:
glm1a=glm(Chao ~ Treatment+Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
# link function made no difference here
par(mfrow=c(2,2))
plot(glm1a)
par(mfrow=c(1,1))
fitgamma <- simulateResiduals(fittedModel = glm1a)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.10183, p-value = 0.702
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.1993, p-value = 0.456

summary(glm1a)
# glm(formula = Chao ~ Treatment + Year, family = Gamma(link = "log"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.89054  -0.32263  -0.07817   0.15871   1.02276  
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          3.58026    0.16247  22.036   <2e-16 ***
# TreatmentHeat        0.02373    0.20551   0.115    0.909    
# TreatmentHeat+Water -0.01210    0.20551  -0.059    0.953    
# TreatmentWater       0.15664    0.20551   0.762    0.450    
# Year2                0.01369    0.14532   0.094    0.925    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.2534127)
# 
# Null deviance: 10.0781  on 47  degrees of freedom
# Residual deviance:  9.8526  on 43  degrees of freedom
# AIC: 412.4
# 
# Number of Fisher Scoring iterations: 6

anova(glm1a, test="F")
# Analysis of Deviance Table 
# Model: Gamma, link: log 
# Response: Chao 
# Terms added sequentially (first to last) 
# 
#           Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                         47     9.9317              
# Treatment  3 0.223309        44     9.8548 0.2937 0.8297
# Year       1 0.002222        43     9.8526 0.0088 0.9258

drop1(glm1a, test="F")
# Single term deletions 
# Model:
#   Chao ~ Treatment + Year
#           Df Deviance    AIC F value Pr(>F)
# <none>         9.6830 411.88               
# Treatment  3  10.0756 407.28  0.3244 0.8077
# Year       1   9.8548 410.41  0.0097 0.9220

emmeans(glm1a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean    SE  df asymp.LCL asymp.UCL
# Control      3.59 0.145 Inf      3.30      3.87
# Heat         3.61 0.145 Inf      3.33      3.90
# Heat+Water   3.58 0.145 Inf      3.29      3.86
# Water        3.74 0.145 Inf      3.46      4.03
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE  df z.ratio p.value
# Control - Heat          -0.0237 0.206 Inf -0.115  0.9994 
# Control - (Heat+Water)   0.0121 0.206 Inf  0.059  0.9999 
# Control - Water         -0.1566 0.206 Inf -0.762  0.8714 
# Heat - (Heat+Water)      0.0358 0.206 Inf  0.174  0.9981 
# Heat - Water            -0.1329 0.206 Inf -0.647  0.9167 
# (Heat+Water) - Water    -0.1687 0.206 Inf -0.821  0.8445 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$Chao~Flowers$Treatment,xlab="Treatment",ylab="Visitor richness",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$Chao,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$Chao,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
# Control       Heat Heat+Water      Water 
# 36.12846   36.95536   35.72713   42.26943 
prse
# Control       Heat Heat+Water      Water 
# 4.500832   4.824713   6.174074   6.087314 
prbarplot=barplot(prmean, ylim=c(0,60), xlab="Treatment",ylab="Mean Visitor richness",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(Chao~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Visitor richness", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$Chao, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$Chao,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
# Control     Heat Heat+Water    Water
# 1 36.10556 42.88690   30.85306 40.13469
# 2 36.15136 31.02381   40.60119 44.40417
prse
# Control     Heat Heat+Water     Water
# 1 3.509941 6.630815   5.343891  3.652113
# 2 8.764305 6.661946  11.387039 12.160714
prbarplot=barplot(prmean, beside=T, ylim=c(0,60), xlab="Treatment",ylab="Mean Visitor richness", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 60, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))


tapply(Flowers$Chao, Flowers$Year, na.rm=TRUE, mean)
# 1        2 
# 37.49505 38.04513 
tapply(Flowers$Chao, Flowers$Year, na.rm=TRUE, std.error)
# 1        2 
# 2.495646 4.770342 


#
########################## 5. Visits Per Flower (all species) ####

# quick look at the distribution:
plot(density(Flowers$VisitsOverFlowers, na.rm=TRUE))

lm1=lm(VisitsOverFlowers ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.93888, p-value = 0.0147
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)

# the residuals look dodgy

# try gamma:
glm1=glm(VisitsOverFlowers ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.098261, p-value = 0.706
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.1375, p-value = 0.312

# quick look at inverse gaussian:
glm2=glm(VisitsOverFlowers ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.10562, p-value = 0.6199
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.0203, p-value = 0.848

# gamma and inverse gaussian look similar (both better than gaussian), gamma is probably slightly better

summary(glm1)
# glm(formula = VisitsOverFlowers ~ Treatment * Year, family = Gamma(link = "log"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.83532  -0.16761  -0.06352   0.14194   0.64725  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -2.83806    0.13030 -21.780   <2e-16 ***
# TreatmentHeat              0.22914    0.18428   1.243   0.2210    
# TreatmentHeat+Water        0.33345    0.18428   1.809   0.0779 .  
# TreatmentWater             0.07782    0.18428   0.422   0.6751    
# Year2                     -0.11650    0.18428  -0.632   0.5309    
# TreatmentHeat:Year2        0.26352    0.26061   1.011   0.3180    
# TreatmentHeat+Water:Year2  0.22311    0.26061   0.856   0.3970    
# TreatmentWater:Year2      -0.22841    0.26061  -0.876   0.3860    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.1018752)
# 
# Null deviance: 6.7175  on 47  degrees of freedom
# Residual deviance: 4.1190  on 40  degrees of freedom
# AIC: -226.66
# 
# Number of Fisher Scoring iterations: 4

anova(glm1, test="F")
# Analysis of Deviance Table 
# Model: Gamma, link: log 
# Response: VisitsOverFlowers 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                              47     6.7175                     
# Treatment       3  2.10375        44     4.6137 6.8834 0.0007593 ***
# Year            1  0.03172        43     4.5820 0.3114 0.5799627    
# Treatment:Year  3  0.46298        40     4.1190 1.5149 0.2253880    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# the interaction isn't significant, will re-run without it:
glm1a=glm(VisitsOverFlowers ~ Treatment+Year, data=Flowers, family=Gamma(link="identity"), na.action = "na.omit")
# link function affected the residuals here, identity gave the best fit
par(mfrow=c(2,2))
plot(glm1a)
par(mfrow=c(1,1))
fitgamma <- simulateResiduals(fittedModel = glm1a)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.094059, p-value = 0.7541
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.1316, p-value = 0.312

summary(glm1a)
# glm(formula = VisitsOverFlowers ~ Treatment + Year, family = Gamma(link = "identity"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.90347  -0.21110  -0.02552   0.13777   0.75793  
# 
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.058149   0.006182   9.406 5.36e-12 ***
# TreatmentHeat        0.024781   0.009157   2.706  0.00972 ** 
# TreatmentHeat+Water  0.031512   0.009687   3.253  0.00223 ** 
# TreatmentWater      -0.001918   0.007220  -0.266  0.79181    
# Year2               -0.005730   0.006075  -0.943  0.35087    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.1068213)
# 
# Null deviance: 6.7175  on 47  degrees of freedom
# Residual deviance: 4.5214  on 43  degrees of freedom
# AIC: -228.12
# 
# Number of Fisher Scoring iterations: 7

anova(glm1a, test="F")
# Analysis of Deviance Table 
# Model: Gamma, link: identity 
# Response: VisitsOverFlowers 
# Terms added sequentially (first to last) 
# 
#           Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                         47     6.7175                     
# Treatment  3   2.1037        44     4.6137 6.5647 0.0009411 ***
# Year       1   0.0923        43     4.5214 0.8640 0.3578013    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm1a, test="F")
# Single term deletions
# 
# Model:
#   VisitsOverFlowers ~ Treatment + Year
#           Df Deviance     AIC F value    Pr(>F)    
# <none>         4.5214 -228.12                      
# Treatment  3   6.7149 -213.59  6.9535 0.0006441 ***
# Year       1   4.6137 -229.25  0.8778 0.3540441    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm1a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean      SE  df asymp.LCL asymp.UCL
# Control    0.0553 0.00520 Inf    0.0423    0.0682
# Heat       0.0801 0.00754 Inf    0.0613    0.0989
# Heat+Water 0.0868 0.00818 Inf    0.0664    0.1072
# Water      0.0534 0.00502 Inf    0.0409    0.0659
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate      SE  df z.ratio p.value
# Control - Heat       -0.02478 0.00916 Inf -2.706  0.0344 
# Control - Heat+Water -0.03151 0.00969 Inf -3.253  0.0063 
# Control - Water       0.00192 0.00722 Inf  0.266  0.9934 
# Heat - Heat+Water    -0.00673 0.01112 Inf -0.605  0.9305 
# Heat - Water          0.02670 0.00905 Inf  2.948  0.0168 
# Heat+Water - Water    0.03343 0.00959 Inf  3.485  0.0028 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 4 estimates 

# quick graphs of raw data:
boxplot(Flowers$VisitsOverFlowers~Flowers$Treatment,xlab="Treatment",ylab="Visits per flower",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$VisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$VisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 0.05532061 0.07944343 0.08630385 0.05404759 
prse
#   Control       Heat Heat+Water      Water 
# 0.003978982 0.007685445 0.010218264 0.004584026 
prbarplot=barplot(prmean, ylim=c(0,0.12), xlab="Treatment",ylab="Mean Visits per flower",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(VisitsOverFlowers~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Visits per flower", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$VisitsOverFlowers, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$VisitsOverFlowers,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 0.05853934 0.07361396 0.08170779 0.06327707
# 2 0.05210189 0.08527290 0.09089990 0.04481811
prse
#    Control     Heat Heat+Water    Water 
# 1 0.005934907 0.005550228 0.01724219 0.007272230
# 2 0.005504080 0.014679666 0.01239655 0.002345174
prbarplot=barplot(prmean, beside=T, ylim=c(0,0.12), xlab="Treatment",ylab="Mean Visits per flower", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1,0.12, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



########################## 6. Visits Per C. cyanus Flower ####

# quick look at the distribution:
plot(density(Flowers$CFVisitsOverFlowers, na.rm=TRUE))

lm1=lm(CFVisitsOverFlowers ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.96854, p-value = 0.2684
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)

# these residuals don't look great, but the data contains true 0s, so can't use Gamma

# transformations all made the residuals worse

# gaussian is the best 

summary(lm1)
# lm(formula = CFVisitsOverFlowers ~ Treatment * Year, data = Flowers)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.277089 -0.039322 -0.008641  0.060182  0.222911 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                0.12154    0.04774   2.546   0.0153 *
# TreatmentHeat              0.04680    0.06752   0.693   0.4927  
# TreatmentHeat+Water        0.10563    0.06752   1.564   0.1265  
# TreatmentWater             0.01681    0.06752   0.249   0.8048  
# Year2                      0.15555    0.07082   2.197   0.0346 *
# TreatmentHeat:Year2       -0.16025    0.10015  -1.600   0.1183  
# TreatmentHeat+Water:Year2 -0.14939    0.10015  -1.492   0.1445  
# TreatmentWater:Year2      -0.13452    0.10015  -1.343   0.1876  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1169 on 36 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1815,	Adjusted R-squared:  0.02231 
# F-statistic:  1.14 on 7 and 36 DF,  p-value: 0.3606

anova(lm1, test="F")
# Analysis of Variance Table 
# Response: CFVisitsOverFlowers
#                Df  Sum Sq  Mean Sq F value Pr(>F)
# Treatment       3 0.04180 0.013934  1.0188 0.3957
# Year            1 0.02161 0.021610  1.5801 0.2168
# Treatment:Year  3 0.04575 0.015249  1.1150 0.3558
# Residuals      36 0.49236 0.013677      

# the interaction is not significant, will re-run without it:
lm1a=lm(CFVisitsOverFlowers ~ Treatment+Year, data=Flowers)
shapiro.test(residuals(lm1a))
# data:  residuals(lm1)
# W = 0.97999, p-value = 0.6332
par(mfrow=c(2,2))
plot(lm1a)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1a, n = 250)
plot(fitg, asFactor = TRUE)

summary(lm1a)
# lm(formula = CFVisitsOverFlowers ~ Treatment + Year, data = Flowers) 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.254243 -0.063577 -0.002123  0.055130  0.283479 
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.17201    0.03893   4.418 7.69e-05 ***
# TreatmentHeat       -0.02604    0.05009  -0.520    0.606    
# TreatmentHeat+Water  0.03772    0.05009   0.753    0.456    
# TreatmentWater      -0.04434    0.05009  -0.885    0.381    
# Year2                0.04451    0.03556   1.251    0.218    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1175 on 39 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1054,	Adjusted R-squared:  0.01367 
# F-statistic: 1.149 on 4 and 39 DF,  p-value: 0.3481

anova(lm1a, test="F")
# Analysis of Variance Table 
# Response: CFVisitsOverFlowers
# Df  Sum Sq  Mean Sq F value Pr(>F)
# Treatment  3 0.04180 0.013934  1.0099 0.3987
# Year       1 0.02161 0.021610  1.5662 0.2182
# Residuals 39 0.53811 0.013798 

drop1(lm1a, test = "F")
# Single term deletions
# Model:
#   CFVisitsOverFlowers ~ Treatment + Year
#           Df Sum of Sq     RSS     AIC F value Pr(>F)
# <none>                 0.53811 -183.77               
# Treatment  3  0.041803 0.57991 -186.48  1.0099 0.3987
# Year       1  0.021610 0.55972 -184.04  1.5662 0.2182

emmeans(lm1a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE df lower.CL upper.CL
# Control     0.194 0.0355 39   0.1017    0.287
# Heat        0.168 0.0355 39   0.0757    0.261
# Heat+Water  0.232 0.0355 39   0.1394    0.325
# Water       0.150 0.0355 39   0.0574    0.243
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE df t.ratio p.value
# Control - Heat         0.0260 0.0501 39  0.520  0.9538 
# Control - Heat+Water  -0.0377 0.0501 39 -0.753  0.8748 
# Control - Water        0.0443 0.0501 39  0.885  0.8125 
# Heat - Heat+Water     -0.0638 0.0501 39 -1.273  0.5852 
# Heat - Water           0.0183 0.0501 39  0.365  0.9831 
# Heat+Water - Water     0.0821 0.0501 39  1.638  0.3697 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 4 estimates

# quick graphs of raw data:
boxplot(Flowers$CFVisitsOverFlowers~Flowers$Treatment,xlab="Treatment",ylab="Visits per CF flower",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$CFVisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$CFVisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
#  0.1922444  0.1662016  0.2299663  0.1479081 
prse
#   Control       Heat Heat+Water      Water 
# 0.05082905 0.02215448 0.03732374 0.02491716 
prbarplot=barplot(prmean, ylim=c(0,0.3), xlab="Treatment",ylab="Mean Visits per CF flower",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(CFVisitsOverFlowers~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Visits per CF flower", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

Flowers2 = Flowers[-which(is.na(Flowers$CFVisitsOverFlowers)),]
prmean=tapply(Flowers2$CFVisitsOverFlowers,list(Flowers2$Year, Flowers2$Treatment), mean)
prse=tapply(Flowers2$CFVisitsOverFlowers,list(Flowers2$Year, Flowers2$Treatment), std.error)
prmean
#    Control      Heat  Heat+Water    Water 
# 1 0.1215405 0.16834    0.2271659 0.1383505
# 2 0.2770889 0.1636354  0.2333269 0.1593773
prse
#    Control     Heat Heat+Water    Water 
# 1 0.02487240 0.01767129 0.04645639 0.02144895
# 2 0.09991218 0.04720259 0.06647834 0.05161765
prbarplot=barplot(prmean, beside=T, ylim=c(0,0.4), xlab="Treatment",ylab="Mean Visits per CF flower", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(9, 0.4, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))

prmean=tapply(Flowers$CFVisitsOverFlowers,Flowers$Year, na.rm=TRUE, mean)
prse=tapply(Flowers$CFVisitsOverFlowers,Flowers$Year, na.rm=TRUE, std.error)
prmean
# 1         2 
# 0.1638492 0.2083571 
prse
# 1          2 
# 0.01622371 0.03381912 


# Worth looking at the number of flowers as well:
tapply(Flowers$CFFlowers, Flowers$Treatment, na.rm=TRUE, mean)
# Control       Heat Heat+Water      Water 
# 124.83333  107.41667   93.58333  106.00000 
tapply(Flowers$CFFlowers, Flowers$Treatment, na.rm=TRUE, std.error)
# Control       Heat Heat+Water      Water 
# 37.31354   18.67970   19.10000   24.97969 



#
########################## 7. Visits Per G. segetum Flower ####

# quick look at the distribution:
plot(density(Flowers$CMVisitsOverFlowers, na.rm=TRUE))

lm1=lm(CMVisitsOverFlowers ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.94001, p-value = 0.01625
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)

# gaussian doesn't look too bad, worth trying gamma:
glm1=glm(CMVisitsOverFlowers ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.13029, p-value = 0.3579
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.0934, p-value = 0.424

# quick look at inverse gaussian:
glm2=glm(CMVisitsOverFlowers ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.10137, p-value = 0.6697
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.90143, p-value = 0.664

# inverse gaussian looks the best fit

summary(glm2)
# glm(formula = CMVisitsOverFlowers ~ Treatment * Year, family = inverse.gaussian(link = "log"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.7915  -0.7987  -0.3226   0.5272   2.0556  
# 
# Coefficients:
#                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               -2.38449    0.12601 -18.923   <2e-16 ***
# TreatmentHeat              0.50186    0.20520   2.446   0.0190 *  
# TreatmentHeat+Water        0.52918    0.20696   2.557   0.0145 *  
# TreatmentWater            -0.08643    0.17448  -0.495   0.6230    
# Year2                      0.03627    0.17984   0.202   0.8412    
# TreatmentHeat:Year2       -0.40978    0.27681  -1.480   0.1466    
# TreatmentHeat+Water:Year2 -0.26061    0.28432  -0.917   0.3648    
# TreatmentWater:Year2       0.11147    0.25257   0.441   0.6613    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 1.034033)
# 
# Null deviance: 59.639  on 47  degrees of freedom
# Residual deviance: 39.085  on 40  degrees of freedom
# AIC: -178.93
# 
# Number of Fisher Scoring iterations: 5

anova(glm2, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: log 
# Response: CMVisitsOverFlowers 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev      F   Pr(>F)   
# NULL                              47     59.639                   
# Treatment       3  15.4198        44     44.219 4.9708 0.005026 **
# Year            1   0.5188        43     43.700 0.5018 0.482836   
# Treatment:Year  3   4.6156        40     39.085 1.4879 0.232418   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# the interaction isn't significant, will re-run without it:
glm2a=glm(CMVisitsOverFlowers ~ Treatment+Year, data=Flowers, family=inverse.gaussian(link="identity"), na.action = "na.omit")
# link function slightly affected the residuals here, identity gave the best fit
par(mfrow=c(2,2))
plot(glm2a)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2a),
                      fitted=fitted.values(glm2a))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2a)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.089873, p-value = 0.8001
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.93373, p-value = 0.768

summary(glm2a)
# glm(formula = CMVisitsOverFlowers ~ Treatment + Year, family = inverse.gaussian(link = "identity"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.7869  -0.8970  -0.2291   0.6121   1.8761  
# 
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.095981   0.010399   9.229 9.29e-12 ***
# TreatmentHeat        0.033429   0.016388   2.040   0.0475 *  
# TreatmentHeat+Water  0.046063   0.018206   2.530   0.0151 *  
# TreatmentWater      -0.002283   0.012179  -0.187   0.8522    
# Year2               -0.003845   0.010536  -0.365   0.7170    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 1.111582)
# 
# Null deviance: 59.639  on 47  degrees of freedom
# Residual deviance: 44.083  on 43  degrees of freedom
# AIC: -179.16
# 
# Number of Fisher Scoring iterations: 8

anova(glm2a, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: identity 
# Response: CMVisitsOverFlowers 
# Terms added sequentially (first to last) 
# 
#           Df Deviance Resid. Df Resid. Dev      F   Pr(>F)   
# NULL                         47     59.639                   
# Treatment  3  15.4198        44     44.219 4.6240 0.006861 **
# Year       1   0.1358        43     44.083 0.1221 0.728440   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm2a, test="F")
# Single term deletions 
# Model:
#   CMVisitsOverFlowers ~ Treatment + Year
#           Df Deviance     AIC F value   Pr(>F)   
# <none>         44.083 -179.16                    
# Treatment  3   57.661 -172.94  4.4145 0.008583 **
# Year       1   44.219 -181.03  0.1324 0.717716   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm2a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean      SE  df asymp.LCL asymp.UCL
# Control    0.0941 0.00877 Inf    0.0722     0.116
# Heat       0.1275 0.01385 Inf    0.0930     0.162
# Heat+Water 0.1401 0.01596 Inf    0.1004     0.180
# Water      0.0918 0.00846 Inf    0.0707     0.113
# 
# Results are averaged over the levels of: Year 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat       -0.03343 0.0164 Inf -2.040  0.1734 
# Control - Heat+Water -0.04606 0.0182 Inf -2.530  0.0554 
# Control - Water       0.00228 0.0122 Inf  0.187  0.9977 
# Heat - Heat+Water    -0.01263 0.0211 Inf -0.598  0.9327 
# Heat - Water          0.03571 0.0162 Inf  2.202  0.1227 
# Heat+Water - Water    0.04835 0.0181 Inf  2.678  0.0372 
# 
# Results are averaged over the levels of: Year 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$CMVisitsOverFlowers~Flowers$Treatment,xlab="Treatment",ylab="Visits per CM flower",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$CMVisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$CMVisitsOverFlowers,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#    Control       Heat Heat+Water      Water 
# 0.09383724 0.12847203 0.14068908 0.09123388 
prse
#     Control        Heat  Heat+Water      Water 
# 0.010222762 0.012853468 0.015792312 0.007315781 
prbarplot=barplot(prmean, ylim=c(0,0.2), xlab="Treatment",ylab="Mean Visits per CM flower",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(CMVisitsOverFlowers~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Visits per CM flower", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$CMVisitsOverFlowers, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$CMVisitsOverFlowers,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#      Control      Heat Heat+Water      Water
# 1 0.09213557 0.1521893  0.1564041 0.08450645
# 2 0.09553890 0.1047548  0.1249741 0.09796132
prse
#      Control       Heat Heat+Water      Water
# 1 0.01786517 0.01778096 0.02815778 0.00511507
# 2 0.01181097 0.01362955 0.01434234 0.01382836
prbarplot=barplot(prmean, beside=T, ylim=c(0,0.2), xlab="Treatment",ylab="Mean Visits per CM flower", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 0.2, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))

# Number of flowers:
tapply(Flowers$CMFlowers,Flowers$Treatment, na.rm=TRUE, mean)
# Control       Heat Heat+Water      Water 
# 706.8333   527.0000   437.2500   640.0833 
tapply(Flowers$CMFlowers,Flowers$Treatment, na.rm=TRUE, std.error)
# Control       Heat Heat+Water      Water 
# 90.41578  108.41789  128.03469   78.28943 



#
########################## 8. Visits Per L. purpureum Flower ####

# this is too data sparse to model




########################## 9. Visits Per V. persica Flower ####

# this is too data sparse to model




########################## 10. Visits Per S. media Flower ####

# this is too data sparse to model




########################## 11. Diet Breadth ####

# quick look at the distribution:
plot(density(Flowers$MeanDiet, na.rm=TRUE))

lm1=lm(MeanDiet ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.95726, p-value = 0.07827
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)

# gaussian doesn't look too bad, worth a quick look at other families

# gamma:
glm1=glm(MeanDiet ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.099507, p-value = 0.6915
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.0264, p-value = 0.784

# quick look at inverse gaussian:
glm2=glm(MeanDiet ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.082949, p-value = 0.8688
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.91308, p-value = 0.432

# these are all fairly similar, but the inverse gaussian is the best fit

summary(glm2)
# glm(formula = MeanDiet ~ Treatment * Year, family = inverse.gaussian(link = "log"), 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.16052  -0.07389  -0.01749   0.07415   0.19063  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.32954    0.05005   6.584 7.14e-08 ***
# TreatmentHeat             -0.05454    0.06983  -0.781   0.4394    
# TreatmentHeat+Water       -0.02090    0.07041  -0.297   0.7681    
# TreatmentWater             0.03762    0.07145   0.526   0.6015    
# Year2                     -0.13359    0.06853  -1.949   0.0583 .  
# TreatmentHeat:Year2       -0.01136    0.09550  -0.119   0.9059    
# TreatmentHeat+Water:Year2  0.05096    0.09699   0.525   0.6022    
# TreatmentWater:Year2      -0.01553    0.09766  -0.159   0.8745    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 0.01080961)
# 
# Null deviance: 0.61972  on 47  degrees of freedom
# Residual deviance: 0.42150  on 40  degrees of freedom
# AIC: -36.974
# 
# Number of Fisher Scoring iterations: 4

anova(glm2, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: log 
# Response: MeanDiet 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                              47    0.61972                      
# Treatment       3 0.040228        44    0.57949  1.2405 0.3077595    
# Year            1 0.151514        43    0.42798 14.0166 0.0005702 ***
# Treatment:Year  3 0.006480        40    0.42150  0.1998 0.8958794    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# the interaction isn't significant, will re-run without it:
glm2a=glm(MeanDiet ~ Treatment+Year, data=Flowers, family=inverse.gaussian, na.action = "na.omit")
sapply(c("log", "1/mu^2", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glm2a, family=inverse.gaussian(link = x)))
  m_AIC = AIC(update(glm2a,  family=inverse.gaussian(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#            log   1/mu^2   inverse  identity
# [1,]  27.12105  27.0447  27.08854  27.13945
# [2,] -42.24210 -42.0894 -42.17707 -42.27889
# link function makes very little difference here. Residual plots slightly better for the canonical (1/mu^2)
par(mfrow=c(2,2))
plot(glm2a)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2a),
                      fitted=fitted.values(glm2a))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2a)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.084861, p-value = 0.8509
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.94909, p-value = 0.728

summary(glm2a)
# glm(formula = MeanDiet ~ Treatment + Year, family = inverse.gaussian, 
#     data = Flowers, na.action = "na.omit")
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.16214  -0.07618  -0.01557   0.07134   0.18484  
# 
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.519762   0.042138  12.335 1.03e-15 ***
# TreatmentHeat        0.073998   0.057717   1.282 0.206687    
# TreatmentHeat+Water -0.003711   0.054900  -0.068 0.946416    
# TreatmentWater      -0.034275   0.053802  -0.637 0.527462    
# Year2                0.152404   0.039934   3.816 0.000429 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 0.01018414)
# 
# Null deviance: 0.61972  on 47  degrees of freedom
# Residual deviance: 0.42934  on 43  degrees of freedom
# AIC: -42.089
# 
# Number of Fisher Scoring iterations: 4

anova(glm2a, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: 1/mu^2 
# Response: MeanDiet 
# Terms added sequentially (first to last) 
# 
#           Df Deviance Resid. Df Resid. Dev       F    Pr(>F)    
# NULL                         47    0.61972                      
# Treatment  3 0.040228        44    0.57949  1.3167 0.2813130    
# Year       1 0.150150        43    0.42934 14.7435 0.0003996 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm2a, test="F")
# Single term deletions# 
# Model:
#   MeanDiet ~ Treatment + Year
#           Df Deviance     AIC F value    Pr(>F)    
# <none>        0.42934 -42.089                      
# Treatment  3  0.46910 -44.185  1.3274 0.2778987    
# Year       1  0.57949 -29.346 15.0380 0.0003559 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm2a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE  df asymp.LCL asymp.UCL
# Control     0.596 0.0391 Inf     0.499     0.693
# Heat        0.670 0.0428 Inf     0.563     0.776
# Heat+Water  0.592 0.0389 Inf     0.495     0.689
# Water       0.562 0.0374 Inf     0.469     0.655
# 
# Results are averaged over the levels of: Year 
# Results are given on the 1/mu^2 (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat       -0.07400 0.0577 Inf -1.282  0.5744 
# Control - Heat+Water  0.00371 0.0549 Inf  0.068  0.9999 
# Control - Water       0.03428 0.0538 Inf  0.637  0.9200 
# Heat - Heat+Water     0.07771 0.0576 Inf  1.349  0.5314 
# Heat - Water          0.10827 0.0565 Inf  1.915  0.2217 
# Heat+Water - Water    0.03056 0.0537 Inf  0.570  0.9412 
# 
# Results are averaged over the levels of: Year 
# Note: contrasts are still on the 1/mu^2 scale 
# P value adjustment: tukey method for comparing a family of 4 estimates


# quick graphs of raw data:
boxplot(Flowers$MeanDiet~Flowers$Treatment,xlab="Treatment",ylab="Diet breadth",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$MeanDiet,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$MeanDiet,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
#  1.303397   1.227713   1.307577   1.343631  
prse
#   Control       Heat Heat+Water      Water 
# 0.04711339 0.04518152 0.05291375 0.04975135 
prbarplot=barplot(prmean, ylim=c(0,1.6), xlab="Treatment",ylab="Mean Diet breadth",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(MeanDiet~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Diet breadth", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$MeanDiet, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$MeanDiet,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
# Control     Heat Heat+Water    Water
# 1 1.390331 1.316538   1.361574 1.443630
# 2 1.216462 1.138889   1.253579 1.243632
prse
# Control       Heat Heat+Water      Water
# 1 0.06597808 0.05659208 0.08222227 0.03306568
# 2 0.04889201 0.05121969 0.06627679 0.07614197
prbarplot=barplot(prmean, beside=T, ylim=c(0,2), xlab="Treatment",ylab="Mean Diet breadth", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 2, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



########################## 12. Weighted Connectance ####

# quick look at the distribution:
plot(density(Flowers$WConnectance, na.rm=TRUE))

# This is a network descriptor, it is bound by 0 and 1 but cannot be examined with a binomial model as it is not
# a proportion of successes. Therefore, the most appropriate analysis method is to use the Beta family in a regression.

glm1 = glmmTMB(WConnectance ~ Treatment*Year, family=beta_family, data=Flowers)
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitbeta <- simulateResiduals(fittedModel = glm1)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.085667, p-value = 0.8727
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 1.0606, p-value = 0.736

summary(glm1)
# Family: beta  ( logit )
# Formula:          WConnectance ~ Treatment * Year
# Data: Flowers
# 
#    AIC      BIC   logLik deviance df.resid 
# -175.5   -158.6     96.7   -193.5       39 
# 
# 
# Overdispersion parameter for beta family ():  143 
# 
# Conditional model:
#                           Estimate Std. Error z value Pr(>|z|)    
# (Intercept)               -1.81513    0.09730 -18.654  < 2e-16 ***
# TreatmentHeat              0.29362    0.13125   2.237  0.02528 *  
# TreatmentHeat+Water        0.31287    0.13090   2.390  0.01684 *  
# TreatmentWater             0.21159    0.13283   1.593  0.11118    
# Year2                      0.40452    0.12930   3.128  0.00176 ** 
# TreatmentHeat:Year2       -0.10841    0.17618  -0.615  0.53834    
# TreatmentHeat+Water:Year2 -0.02685    0.17496  -0.153  0.87802    
# TreatmentWater:Year2      -0.28341    0.18028  -1.572  0.11594    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


drop1(glm1, test = "Chisq")
# Single term deletions
# 
# Model:
#   WConnectance ~ Treatment * Year
#                Df     AIC    LRT Pr(>Chi)
# <none>            -175.47                
# Treatment:Year  3 -178.46 3.0106     0.39


glmmTMB:::Anova.glmmTMB(glm1)
# Analysis of Deviance Table (Type II Wald chisquare tests) 
# Response: WConnectance
# Chisq Df Pr(>Chisq)    
#   Treatment      15.6867  3   0.001315 ** 
#   Year           23.9826  1  9.721e-07 ***
#   Treatment:Year  3.1159  3   0.374103    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
glmmTMB:::Anova.glmmTMB(glm1, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: WConnectance
#                     Chisq Df Pr(>Chisq)    
#   (Intercept)    347.9887  1  < 2.2e-16 ***
#   Treatment        6.9338  3   0.074039 .  
#   Year             9.7874  1   0.001757 ** 
#   Treatment:Year   3.1159  3   0.374103    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# double check the significance of the interaction term using likelihood ratio test:
glm2 = glmmTMB(WConnectance ~ Treatment + Year, family=beta_family, data=Flowers)
lrtest(glm1, glm2)
# Likelihood ratio test 
# Model 1: WConnectance ~ Treatment * Year
# Model 2: WConnectance ~ Treatment + Year
# #  Df LogLik Df  Chisq Pr(>Chisq)
# 1   9 96.734                     
# 2   6 95.228 -3 3.0106       0.39

# this shows the interaction isn't significant


# quick check of other link functions:
sapply(c("logit", "probit", "cloglog"), function(x) {
  bm_loglik = logLik(update(glm2, family = beta_family(link = x)))
  bm_AIC = AIC(update(glm2, family = beta_family(link = x)))
  df1 = rbind(bm_loglik, bm_AIC)
  df1
})
#           logit     probit    cloglog
# [1,]   95.22832   95.19622   95.25234
# [2,] -178.45664 -178.39244 -178.50468

# these are all basically the same, so I will stick with the canonical (logit).


fitbeta <- simulateResiduals(fittedModel = glm2)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.085833, p-value = 0.8713
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 1.0477, p-value = 0.744

glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)

summary(glm2)
# Family: beta  ( logit )
# Formula:          WConnectance ~ Treatment + Year
# Data: Flowers
# 
# AIC      BIC   logLik deviance df.resid 
# -178.5   -167.2     95.2   -190.5       42 
# 
# 
# Overdispersion parameter for beta family ():  135 
# 
# Conditional model:
#                       Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)         -1.75752    0.07417 -23.696  < 2e-16 ***
#   TreatmentHeat        0.23395    0.09017   2.595 0.009469 ** 
#   TreatmentHeat+Water  0.29736    0.08938   3.327 0.000878 ***
#   TreatmentWater       0.05618    0.09262   0.607 0.544166    
#   Year2                0.30145    0.06337   4.757 1.96e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm2, test = "Chisq")
# Single term deletions# 
# Model:
#   WConnectance ~ Treatment + Year
#           Df     AIC    LRT  Pr(>Chi)    
# <none>       -178.46                     
# Treatment  3 -171.34 13.118  0.004388 ** 
# Year       1 -161.83 18.625 1.592e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

glmmTMB:::Anova.glmmTMB(glm2)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: WConnectance
#              Chisq Df Pr(>Chisq)    
#   Treatment 15.049  3   0.001775 ** 
#   Year      22.632  1  1.962e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


emmeans(glm2, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE df lower.CL upper.CL
# Control     -1.61 0.0661 42    -1.74    -1.47
# Heat        -1.37 0.0616 42    -1.50    -1.25
# Heat+Water  -1.31 0.0604 42    -1.43    -1.19
# Water       -1.55 0.0653 42    -1.68    -1.42
# 
# Results are averaged over the levels of: Year 
# Results are given on the logit (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate     SE df t.ratio p.value
# Control - Heat          -0.2339 0.0902 42 -2.595  0.0602 
# Control - (Heat+Water)  -0.2974 0.0894 42 -3.327  0.0095 
# Control - Water         -0.0562 0.0926 42 -0.607  0.9294 
# Heat - (Heat+Water)     -0.0634 0.0861 42 -0.736  0.8820 
# Heat - Water             0.1778 0.0895 42  1.987  0.2092 
# (Heat+Water) - Water     0.2412 0.0887 42  2.719  0.0451 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$WConnectance~Flowers$Treatment,xlab="Treatment",ylab="Weighted Connectance",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$WConnectance,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$WConnectance,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#   Control       Heat Heat+Water      Water 
# 0.1675333  0.2052239  0.2130770  0.1750102 
prse
#     Control        Heat  Heat+Water      Water 
# 0.012101149 0.016610686 0.012450964 0.006994033 
prbarplot=barplot(prmean, ylim=c(0,0.3), xlab="Treatment",ylab="Mean Weighted Connectance",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(WConnectance~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Weighted Connectance", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$WConnectance, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$WConnectance,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#     Control      Heat Heat+Water    Water 
# 1 0.1388407 0.1812097  0.1820920 0.1656831
# 2 0.1962259 0.2292381  0.2440621 0.1843373
prse
#       Control       Heat  Heat+Water      Water 
# 1 0.008919596 0.02035598 0.014094288 0.006836419
# 2 0.015344738 0.02385345 0.009971867 0.011562793
prbarplot=barplot(prmean, beside=T, ylim=c(0,0.3), xlab="Treatment",ylab="Mean Weighted Connectance", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 0.3, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))


###

# I belatedly realised that glmmTMB supports the beta family via the betareg package, which
# means I can validate the model using DHARMa.

# I will do a quick check of the model fit:

glm1 = glmmTMB(WConnectance ~ Treatment*Year, family=beta_family, data=Flowers)
fitbeta <- simulateResiduals(fittedModel = glm1)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.085667, p-value = 0.8727
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 1.0606, p-value = 0.736

# this all looks good.

glm2 = glmmTMB(WConnectance ~ Treatment + Year, family=beta_family, data=Flowers)
fitbeta <- simulateResiduals(fittedModel = glm2)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.085833, p-value = 0.8713
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 1.0477, p-value = 0.744

# it also looks good without the interaction term.




#
########################## 15. Interaction Evenness ####

# quick look at the distribution:
plot(density(Flowers$InteractionEvenness, na.rm=TRUE))

# This is a network descriptor, it is bound by 0 and 1 but cannot be examined with a binomial model as it is not
# a proportion of successes. Therefore, the most appropriate analysis method is to use Beta regression.

glm1 = glmmTMB(InteractionEvenness ~ Treatment*Year, family=beta_family, data=Flowers)
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitbeta <- simulateResiduals(fittedModel = glm1)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.10933, p-value = 0.6146
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 0.96378, p-value = 0.888

summary(glm1)
# Family: beta  ( logit )
# Formula:          InteractionEvenness ~ Treatment * Year
# Data: Flowers
# AIC      BIC   logLik deviance df.resid 
# -121.2   -104.4     69.6   -139.2       39 
# 
# Overdispersion parameter for beta family (): 71.3 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                0.45279    0.09844   4.600 4.23e-06 ***
#   TreatmentHeat              0.25040    0.14163   1.768   0.0771 .  
# TreatmentHeat+Water        0.05657    0.13965   0.405   0.6854    
# TreatmentWater             0.08716    0.13992   0.623   0.5333    
# Year2                     -0.20040    0.13803  -1.452   0.1465    
# TreatmentHeat:Year2        0.08590    0.19860   0.433   0.6654    
# TreatmentHeat+Water:Year2  0.38059    0.19797   1.922   0.0545 .  
# TreatmentWater:Year2      -0.03056    0.19590  -0.156   0.8760    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm1, test = "Chisq")
# Single term deletions
# 
# Model:
#   WConnectance ~ Treatment * Year
#                Df     AIC    LRT Pr(>Chi)
# <none>            -121.24               
# Treatment:Year  3 -122.19 5.049   0.1682


glmmTMB:::Anova.glmmTMB(glm1)
# Analysis of Deviance Table (Type II Wald chisquare tests) 
# Response: WConnectance
# Chisq Df Pr(>Chisq)    
# Treatment      11.8473  3   0.007925 **
# Year            1.8046  1   0.179156   
# Treatment:Year  5.3182  3   0.149923 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
glmmTMB:::Anova.glmmTMB(glm1, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# Response: InteractionEvenness
# Chisq Df Pr(>Chisq)    
# (Intercept)    21.1588  1  4.227e-06 ***
#   Treatment       3.4162  3     0.3318    
# Year            2.1078  1     0.1465    
# Treatment:Year  5.3182  3     0.1499    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# this shows the interaction isn't significant

# quick check of other link functions:
glm2 = glmmTMB(InteractionEvenness ~ Treatment + Year, family=beta_family, data=Flowers)
sapply(c("logit", "probit", "cloglog"), function(x) {
  bm_loglik = logLik(update(glm2, family = beta_family(link = x)))
  bm_AIC = AIC(update(glm2, family = beta_family(link = x)))
  df1 = rbind(bm_loglik, bm_AIC)
  df1
})
#           logit     probit    cloglog
# [1,]   67.09555   67.08433   67.03893
# [2,] -122.19109 -122.16866 -122.07786

# these are all basically the same, so I will stick with the canonical (logit).


fitbeta <- simulateResiduals(fittedModel = glm2)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.095, p-value = 0.7792
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 0.96735, p-value = 0.896

glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)

summary(glm2)
# Family: beta  ( logit )
# Formula:          InteractionEvenness ~ Treatment + Year
# Data: Flowers# 
# AIC      BIC   logLik deviance df.resid 
# -122.2   -111.0     67.1   -134.2       42  
# 
# Overdispersion parameter for beta family (): 64.2 
# 
# Conditional model:
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          0.39873    0.08141   4.898 9.69e-07 ***
# TreatmentHeat        0.29311    0.10446   2.806  0.00502 ** 
# TreatmentHeat+Water  0.24876    0.10422   2.387  0.01699 *  
# TreatmentWater       0.07154    0.10295   0.695  0.48713    
# Year2               -0.09423    0.07388  -1.275  0.20216    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

drop1(glm2, test = "Chisq")
# Single term deletions# 
# Model:
#   WConnectance ~ Treatment + Year
#           Df     AIC    LRT  Pr(>Chi)    
# <none>       -178.46                     
# Treatment  3 -118.45 9.7432  0.02088 *
# Year       1 -122.59 1.6001  0.20588 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

glmmTMB:::Anova.glmmTMB(glm2)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: WConnectance
#              Chisq Df Pr(>Chisq)    
# Treatment 10.7956  3    0.01288 *
# Year       1.6267  1    0.20216  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


emmeans(glm2, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE df lower.CL upper.CL
# Control     0.352 0.0726 42    0.205    0.498
# Heat        0.645 0.0752 42    0.493    0.796
# Heat+Water  0.600 0.0748 42    0.449    0.751
# Water       0.423 0.0730 42    0.276    0.571
# 
# Results are averaged over the levels of: Year 
# Results are given on the logit (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast               estimate    SE df t.ratio p.value
# Control - Heat          -0.2931 0.104 42 -2.806  0.0366 
# Control - (Heat+Water)  -0.2488 0.104 42 -2.387  0.0953 
# Control - Water         -0.0715 0.103 42 -0.695  0.8985 
# Heat - (Heat+Water)      0.0443 0.106 42  0.418  0.9750 
# Heat - Water             0.2216 0.105 42  2.114  0.1652 
# (Heat+Water) - Water     0.1772 0.105 42  1.695  0.3389 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$InteractionEvenness~Flowers$Treatment,xlab="Treatment",ylab="Interaction Evenness",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$InteractionEvenness,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$InteractionEvenness,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
#  0.5869513  0.6530873  0.6458516  0.6054604  
prse
#   Control       Heat Heat+Water      Water 
# 0.02032762 0.02338910 0.01478724 0.01095747 
prbarplot=barplot(prmean, ylim=c(0,0.8), xlab="Treatment",ylab="Mean Interaction Evenness",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(InteractionEvenness~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Interaction Evenness", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$InteractionEvenness, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$InteractionEvenness,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 0.6120339 0.6679866  0.6260683 0.6333668
# 2 0.5618687 0.6381881  0.6656349 0.5775539
prse
#    Control     Heat Heat+Water    Water 
# 1 0.01822185 0.02812902 0.01143049 0.01066503
# 2 0.03513463 0.03907664 0.02597896 0.01015088
prbarplot=barplot(prmean, beside=T, ylim=c(0,0.8), xlab="Treatment",ylab="Mean Interaction Evenness", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 0.8, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))


###

# I belatedly realised that glmmTMB supports the beta family via the betareg package, which
# means I can validate the model using DHARMa.

# I will do a quick check of the model fit:

glm1 = glmmTMB(InteractionEvenness ~ Treatment*Year, family=beta_family, data=Flowers)
fitbeta <- simulateResiduals(fittedModel = glm1)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.10933, p-value = 0.6146
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 0.96378, p-value = 0.888

# this all looks fine.

glm2 = glmmTMB(InteractionEvenness ~ Treatment + Year, family=beta_family, data=Flowers)
fitbeta <- simulateResiduals(fittedModel = glm2)
plot(fitbeta, asFactor = TRUE)
testUniformity(fitbeta)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.095, p-value = 0.7792
testDispersion(fitbeta)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# dispersion = 0.96735, p-value = 0.896

# it also looks fine without the interaction term.




#
########################## 14. Generality ####

# quick look at the distribution:
plot(density(Flowers$GeneralityHL, na.rm=TRUE))

# this is a network descriptor, it is bound by 0 and continuous. A linear model could provide a good fit
# for the data

lm1=lm(GeneralityHL ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.96385, p-value = 0.1444
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)
lm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(lm1),
                      fitted=fitted.values(lm1))
ggplot(lm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(lm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)


# this isn't too bad, but worth trying other families 

# try gamma:
glm1=glm(GeneralityHL ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.084739, p-value = 0.8521
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.0106, p-value = 0.92

# quick look at inverse gaussian:
glm2=glm(GeneralityHL ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.061493, p-value = 0.9883
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.89235, p-value = 0.36

# the Gamma and inverse gaussian are both better than gaussian, and the inverse gaussian looks slightly better 
# than gamma

summary(glm2)
# glm(formula = GeneralityHL ~ Treatment * Year, family = inverse.gaussian(link = "log"), 
#     data = Flowers, na.action = "na.omit") 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.24581  -0.10196  -0.02162   0.08992   0.22308  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                0.42561    0.07120   5.977 5.08e-07 ***
#   TreatmentHeat              0.02185    0.10125   0.216   0.8302    
# TreatmentHeat+Water       -0.01526    0.10031  -0.152   0.8798    
# TreatmentWater             0.05006    0.10198   0.491   0.6262    
# Year2                     -0.16978    0.09668  -1.756   0.0867 .  
# TreatmentHeat:Year2       -0.01207    0.13730  -0.088   0.9304    
# TreatmentHeat+Water:Year2  0.08262    0.13754   0.601   0.5514    
# TreatmentWater:Year2      -0.04630    0.13774  -0.336   0.7385    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 0.01987485)
# 
# Null deviance: 1.0501  on 47  degrees of freedom
# Residual deviance: 0.8027  on 40  degrees of freedom
# AIC: 7.7293
# 
# Number of Fisher Scoring iterations: 4

anova(glm2, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: log 
# Response: GeneralityHL 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev       F   Pr(>F)   
# NULL                              47    1.05007                    
# Treatment       3 0.003991        44    1.04608  0.0669 0.977141   
# Year            1 0.224829        43    0.82125 11.3122 0.001707 **
# Treatment:Year  3 0.018556        40    0.80270  0.3112 0.817145   

# the interaction isn't significant, so will re-run the model without it:
glm2a=glm(GeneralityHL ~ Treatment+Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
sapply(c("log", "1/mu^2", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glm2a, family=inverse.gaussian(link = x)))
  m_AIC = AIC(update(glm2a,  family=inverse.gaussian(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#           log   1/mu^2 inverse identity
# [1,] 4.586859 4.570351 4.57567 4.602905
# [2,] 2.826282 2.859299 2.84866 2.794190
# link function barely affected the model fit or residual plots. As this is a weighted measure that includes
# sample size in the calculation of the value, it is technically an extensive variable. As a result, it is probably best to go
# with the log link.
par(mfrow=c(2,2))
plot(glm2a)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2a),
                      fitted=fitted.values(glm2a))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2a)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.089108, p-value = 0.8082
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.92823, p-value = 0.576

summary(glm2a)
# glm(formula = GeneralityHL ~ Treatment + Year, family = inverse.gaussian(link = "log"), 
#     data = Flowers, na.action = "na.omit") 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.23955  -0.11056  -0.01859   0.09371   0.23197  
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.42248    0.05353   7.892 6.75e-10 ***
# TreatmentHeat        0.01536    0.06666   0.230  0.81885    
# TreatmentHeat+Water  0.03015    0.06691   0.451  0.65454    
# TreatmentWater       0.02535    0.06683   0.379  0.70632    
# Year2               -0.16399    0.04753  -3.450  0.00127 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for inverse.gaussian family taken to be 0.01888626)
# 
# Null deviance: 1.05007  on 47  degrees of freedom
# Residual deviance: 0.82125  on 43  degrees of freedom
# AIC: 2.8263
# 
# Number of Fisher Scoring iterations: 5

anova(glm2a, test="F")
# Analysis of Deviance Table 
# Model: inverse.gaussian, link: log 
# Response: GeneralityHL 
# Terms added sequentially (first to last) 
# 
# Df Deviance Resid. Df Resid. Dev       F   Pr(>F)   
# NULL                         47    1.05007                    
# Treatment  3 0.003991        44    1.04608  0.0704 0.975428   
# Year       1 0.224829        43    0.82125 11.9044 0.001267 **

drop1(glm2a, test="F")
# Single term deletions 
# Model:
#   GeneralityHL ~ Treatment + Year
#           Df Deviance     AIC F value  Pr(>F)   
# <none>        0.82125  2.8263                   
# Treatment  3  0.82573 -2.9366  0.0782 0.97149   
# Year       1  1.04608 12.7307 11.7718 0.00134 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(glm2a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE  df asymp.LCL asymp.UCL
# Control     0.340 0.0470 Inf     0.223     0.458
# Heat        0.356 0.0474 Inf     0.238     0.474
# Heat+Water  0.371 0.0477 Inf     0.252     0.489
# Water       0.366 0.0476 Inf     0.247     0.484
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate     SE  df z.ratio p.value
# Control - Heat       -0.01536 0.0667 Inf -0.230  0.9957 
# Control - Heat+Water -0.03015 0.0669 Inf -0.451  0.9695 
# Control - Water      -0.02535 0.0668 Inf -0.379  0.9814 
# Heat - Heat+Water    -0.01479 0.0672 Inf -0.220  0.9962 
# Heat - Water         -0.00999 0.0671 Inf -0.149  0.9988 
# Heat+Water - Water    0.00480 0.0673 Inf  0.071  0.9999 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 


# quick graphs of raw data:
boxplot(Flowers$GeneralityHL~Flowers$Treatment,xlab="Treatment",ylab="Generality",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$GeneralityHL,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$GeneralityHL,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
# 1.411025   1.434282   1.444424   1.452742 
prse
#   Control       Heat Heat+Water      Water 
# 0.07321947 0.08325854 0.07302512 0.06766948 
prbarplot=barplot(prmean, ylim=c(0,2), xlab="Treatment",ylab="Mean Generality",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(GeneralityHL~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Generality", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$GeneralityHL, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$GeneralityHL,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 1.530522 1.564339   1.507337 1.609090
# 2 1.291528 1.304226   1.381511 1.296395
prse
#    Control     Heat Heat+Water    Water 
# 1 0.09570578 0.11817256 0.09628569 0.04286888
# 2 0.09336701 0.09884369 0.11229243 0.09237284
prbarplot=barplot(prmean, beside=T, ylim=c(0,2), xlab="Treatment",ylab="Mean Generality", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 2, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



########################## 15. Vulnerability ####

# quick look at the distribution:
plot(density(Flowers$VulnerabilityLL, na.rm=TRUE))

# this is a network descriptor, it is bound by 0 and continuous. A linear model could provide a good fit
# for the data

lm1=lm(VulnerabilityLL ~ Treatment*Year, data=Flowers)
shapiro.test(residuals(lm1))
# data:  residuals(lm1)
# W = 0.98437, p-value = 0.7651
par(mfrow=c(2,2))
plot(lm1)
par(mfrow=c(1,1))
fitg <- simulateResiduals(fittedModel = lm1, n = 250)
plot(fitg, asFactor = TRUE)
lm1_df <- data.frame(Treatment=Flowers$Treatment, 
                     Year=Flowers$Year,
                     TrYe=paste0(Flowers$Treatment, Flowers$Year),
                     resid=residuals(lm1),
                     fitted=fitted.values(lm1))
ggplot(lm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(lm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)

# this looks quite good, but worth trying other families that are bound by 0 

# try gamma:
glm1=glm(VulnerabilityLL ~ Treatment*Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm1)
par(mfrow=c(1,1))
glm1_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm1),
                      fitted=fitted.values(glm1))
ggplot(glm1_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm1_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm1)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# D = 0.066298, p-value = 0.975
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 1.0213, p-value = 0.792

# quick look at inverse gaussian:
glm2=glm(VulnerabilityLL ~ Treatment*Year, data=Flowers, family=inverse.gaussian(link="log"), na.action = "na.omit")
par(mfrow=c(2,2))
plot(glm2)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2),
                      fitted=fitted.values(glm2))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitig <- simulateResiduals(fittedModel = glm2)
plot(fitig, asFactor = TRUE)
testUniformity(fitig)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.091925, p-value = 0.7779
testDispersion(fitig)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# ratioObsSim = 0.93165, p-value = 0.616

# they are all pretty good fits, but gamma looks the best

summary(glm1)
# glm(formula = VulnerabilityLL ~ Treatment * Year, family = Gamma(link = "log"), 
#     data = Flowers, na.action = "na.omit") 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4031  -0.1119  -0.0222   0.1016   0.3838  
# 
# Coefficients:
#                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                1.968439   0.083131  23.679   <2e-16 ***
# TreatmentHeat              0.004627   0.117565   0.039    0.969    
# TreatmentHeat+Water       -0.001986   0.117565  -0.017    0.987    
# TreatmentWater             0.048536   0.117565   0.413    0.682    
# Year2                     -0.004099   0.117565  -0.035    0.972    
# TreatmentHeat:Year2        0.098058   0.166262   0.590    0.559    
# TreatmentHeat+Water:Year2  0.116071   0.166262   0.698    0.489    
# TreatmentWater:Year2       0.002666   0.166262   0.016    0.987    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.04146472)
# 
# Null deviance: 1.7841  on 47  degrees of freedom
# Residual deviance: 1.6935  on 40  degrees of freedom
# AIC: 184.88
# 
# Number of Fisher Scoring iterations: 4

anova(glm1, test="F")
# Analysis of Deviance Table 
# Model: Gamma, link: log 
# Response: VulnerabilityLL 
# Terms added sequentially (first to last) 
# 
#                Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                              47     1.7841              
# Treatment       3 0.026456        44     1.7576 0.2127 0.8870
# Year            1 0.030094        43     1.7275 0.7258 0.3993
# Treatment:Year  3 0.034031        40     1.6935 0.2736 0.8441

# interaction not significant, will re-run without
glm2a=glm(VulnerabilityLL ~ Treatment+Year, data=Flowers, family=Gamma(link="log"), na.action = "na.omit")
sapply(c("log", "inverse", "identity"), function(x) {
  m_loglik = logLik(update(glm2a, family=Gamma(link = x)))
  m_AIC = AIC(update(glm2a,  family=Gamma(link = x)))
  df1 = rbind(m_loglik, m_AIC)
})
#            log   inverse  identity
# [1,] -83.91966 -83.90596 -83.93329
# [2,] 179.83933 179.81192 179.86657
# link function barely affected the model fit or residual plots. As this is a weighted measure that includes
# sample size in the calculation of the value, it is technically an extensive variable. As a result, it is probably best to go
# with the log link.
par(mfrow=c(2,2))
plot(glm2a)
par(mfrow=c(1,1))
glm2_df <- data.frame(Treatment=Flowers$Treatment, 
                      Year=Flowers$Year,
                      TrYe=paste0(Flowers$Treatment, Flowers$Year),
                      resid=residuals(glm2a),
                      fitted=fitted.values(glm2a))
ggplot(glm2_df,aes(x=TrYe,y=resid, colour=Treatment))+geom_point()+geom_boxplot(aes(group=TrYe),alpha = 0)
ggplot(glm2_df,aes(x=Treatment,y=resid, colour=Year))+geom_point()+geom_boxplot(aes(group=Treatment),alpha = 0)
fitgamma <- simulateResiduals(fittedModel = glm2a)
plot(fitgamma, asFactor = TRUE)
testUniformity(fitgamma)
# One-sample Kolmogorov-Smirnov test 
# data:  simulationOutput$scaledResiduals
# 0.073504, p-value = 0.9407
testDispersion(fitgamma)
# DHARMa nonparametric dispersion test via sd of residuals fitted vs. simulated
# data:  simulationOutput
# 1.0218, p-value = 0.792

summary(glm2a)
# glm(formula = VulnerabilityLL ~ Treatment + Year, family = Gamma(link = "log"), 
#     data = Flowers, na.action = "na.omit") 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -0.42687  -0.12513  -0.01384   0.10648   0.41843  
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          1.94171    0.06403  30.325   <2e-16 ***
# TreatmentHeat        0.05353    0.08099   0.661    0.512    
# TreatmentHeat+Water  0.05616    0.08099   0.693    0.492    
# TreatmentWater       0.04983    0.08099   0.615    0.542    
# Year2                0.05010    0.05727   0.875    0.387    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 0.03935938)
# 
# Null deviance: 1.7841  on 47  degrees of freedom
# Residual deviance: 1.7275  on 43  degrees of freedom
# AIC: 179.84
# 
# Number of Fisher Scoring iterations: 4

anova(glm2a, test="F")
# Analysis of Deviance Table 
# Model: Gamma, link: log 
# Response: VulnerabilityLL 
# Terms added sequentially (first to last) 
# 
# Df Deviance Resid. Df Resid. Dev      F Pr(>F)
# NULL                         47     1.7841              
# Treatment  3 0.026456        44     1.7576 0.2241 0.8792
# Year       1 0.030094        43     1.7275 0.7646 0.3868

drop1(glm2a, test="F")
# Single term deletions 
# Model:
#   VulnerabilityLL ~ Treatment + Year
#           Df Deviance    AIC F value Pr(>F)
# <none>         1.7275 179.84               
# Treatment  3   1.7530 174.49  0.2112 0.8881
# Year       1   1.7576 178.60  0.7491 0.3916

emmeans(glm2a, pairwise~Treatment, adjust="tukey")
# $emmeans
# Treatment  emmean     SE  df asymp.LCL asymp.UCL
# Control      1.97 0.0573 Inf      1.82      2.11
# Heat         2.02 0.0573 Inf      1.88      2.16
# Heat+Water   2.02 0.0573 Inf      1.88      2.17
# Water        2.02 0.0573 Inf      1.87      2.16
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 4 estimates 
# 
# $contrasts
# contrast             estimate    SE  df z.ratio p.value
# Control - Heat       -0.05353 0.081 Inf -0.661  0.9117 
# Control - Heat+Water -0.05616 0.081 Inf -0.693  0.8997 
# Control - Water      -0.04983 0.081 Inf -0.615  0.9273 
# Heat - Heat+Water    -0.00263 0.081 Inf -0.032  1.0000 
# Heat - Water          0.00370 0.081 Inf  0.046  1.0000 
# Heat+Water - Water    0.00633 0.081 Inf  0.078  0.9998 
# 
# Results are averaged over the levels of: Year 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates


# quick graphs of raw data:
boxplot(Flowers$VulnerabilityLL~Flowers$Treatment,xlab="Treatment",ylab="Vulnerability",col=c("linen","plum","slateblue1","lightblue"))

prmean=tapply(Flowers$VulnerabilityLL,Flowers$Treatment, na.rm=TRUE, mean)
prse=tapply(Flowers$VulnerabilityLL,Flowers$Treatment, na.rm=TRUE, std.error)
prmean
#  Control       Heat Heat+Water      Water 
#  7.144848   7.546990   7.568578   7.510176 
prse
#   Control       Heat Heat+Water      Water 
# 0.4274734  0.4358455  0.4996222  0.3286017 
prbarplot=barplot(prmean, ylim=c(0,10), xlab="Treatment",ylab="Mean Vulnerability",col=c("linen","plum","slateblue1","lightblue"),xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)


boxplot(VulnerabilityLL~Year*Treatment, data=Flowers,xlab="Treatment",ylab="Vulnerability", col=c("grey40","grey90"), boxwex=0.5)
for(i in seq(0.5 , 20 , 2)) {abline(v=i,lty=1, col="grey")}

prmean=tapply(Flowers$VulnerabilityLL, list(Flowers$Year, Flowers$Treatment), mean)
prse=tapply(Flowers$VulnerabilityLL,list(Flowers$Year, Flowers$Treatment), std.error)
prmean
#    Control      Heat Heat+Water    Water 
# 1 7.159493 7.19270   7.145287 7.515560
# 2 7.130203 7.90128   7.991870 7.504791
prse
#    Control     Heat Heat+Water    Water 
# 1 0.4357894 0.5652895  0.5897851 0.5488190
# 2 0.7836001 0.6826924  0.8239064 0.4170062
prbarplot=barplot(prmean, beside=T, ylim=c(0,10), xlab="Treatment",ylab="Mean Vulnerability", col=c("grey40","grey90"), xpd=F)
arrows(prbarplot,prmean+prse,prbarplot,prmean-prse,angle=90,code=3,length=0.1)
abline(h=0)
legend(1, 10, legend=c("Year 1", "Year 2"), fill=c("grey40", "grey90"))



