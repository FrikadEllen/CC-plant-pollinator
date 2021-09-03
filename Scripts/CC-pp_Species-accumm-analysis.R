.libPaths(c("E:\\Rlib", .libPaths()))

library(vegan)



#### Data

raw_int=read.csv("Data/Raw_InteractionData_BothYears.csv", header=TRUE,sep=",", as.is=TRUE, na.strings="(null)")
raw_int$Record_ID = 1:length(raw_int$Specimen_ID)
raw_int_2014 = raw_int[which(raw_int$Year == 2014),]
raw_int_2015 = raw_int[which(raw_int$Year == 2015),]

head(raw_int)
#   Year Plot Round            Plant      Specimen_ID Basic_Type  Treatment Record_ID
# 1 2014    2     1 Veronica_persica Sphaerophoria_sp   Hoverfly Heat+Water         1
# 2 2014    2     1 Lamium_purpureum Bombus_pascuorum  Bumblebee Heat+Water         2
# 3 2014   16     1  Stellaria_media  Eristalis_tenax   Hoverfly    Control         3
# 4 2014    1     1  Stellaria_media  Anthomyiidae_sp        Fly    Control         4
# 5 2014   11     1 Veronica_persica  Eristalis_tenax   Hoverfly       Heat         5
# 6 2014   23     1 Veronica_persica    Tachinidae_sp        Fly    Control         6


####

# 1. Each treatment
# 2. Each plot - 2014
# 3. Each plot - 2015


####

########################## 1. Each treatment ####

#### start with 2014

treat_2014_ls = split(raw_int_2014, raw_int_2014$Treatment)

treat_mat_2014_ls = lapply(treat_2014_ls, function(x) {
  mat1 = tapply(x$Record_ID, list(x$Round, x$Specimen_ID), length)
  mat1[is.na(mat1)] = 0
  mat1
})

head(treat_mat_2014_ls[1])
names(treat_mat_2014_ls)
# [1] "Control"    "Heat"       "Heat+Water" "Water"   

plot(specaccum(treat_mat_2014_ls[[1]], method="random"))
specpool(treat_mat_2014_ls[[1]])

plot(specaccum(treat_mat_2014_ls[[2]], method="random"))
specpool(treat_mat_2014_ls[[2]])

plot(specaccum(treat_mat_2014_ls[[3]], method="random"))
specpool(treat_mat_2014_ls[[3]])

plot(specaccum(treat_mat_2014_ls[[4]], method="random"))
specpool(treat_mat_2014_ls[[4]])
# none reach asymptote.


### 2015

treat_2015_ls = split(raw_int_2015, raw_int_2015$Treatment)

treat_mat_2015_ls = lapply(treat_2015_ls, function(x) {
  mat1 = tapply(x$Record_ID, list(x$Round, x$Specimen_ID), length)
  mat1[is.na(mat1)] = 0
  mat1
})

head(treat_mat_2015_ls[1])
names(treat_mat_2015_ls)
# [1] "Control"    "Heat"       "Heat+Water" "Water"   

plot(specaccum(treat_mat_2015_ls[[1]], method="random"))
specpool(treat_mat_2015_ls[[1]])

plot(specaccum(treat_mat_2015_ls[[2]], method="random"))
specpool(treat_mat_2015_ls[[2]])

plot(specaccum(treat_mat_2015_ls[[3]], method="random"))
specpool(treat_mat_2015_ls[[3]])

plot(specaccum(treat_mat_2015_ls[[4]], method="random"))
specpool(treat_mat_2015_ls[[4]])
# none reach asymptote.



#
########################## 2. Each plot - 2014 ####

raw_int_2014$Round = factor(raw_int_2014$Round)

plot_2014_ls = split(raw_int_2014, raw_int_2014$Plot)

plot_mat_2014_ls = lapply(plot_2014_ls, function(x) {
  mat1 = tapply(x$Record_ID, list(x$Round, x$Specimen_ID), length)
  mat1[is.na(mat1)] = 0
  mat1
})

head(plot_mat_2014_ls[1])
names(plot_mat_2014_ls)
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24"

plot(specaccum(plot_mat_2014_ls[[1]], method="random"))
# cycling through all plots shows none reach asymptote.


# calculate extrapolated richness values:
plot_2014_ex = Reduce(rbind,
                      Map(function(x, y) {
                        df1 = specpool(x)
                        df1$Plot = paste0(y)
                        df1$Year = "2014"
                        rownames(df1) = NULL
                        df1[, c(11, 10, 1:9)]
                        }, plot_mat_2014_ls, names(plot_mat_2014_ls)))
plot_2014_ex
#    Year Plot Species     chao   chao.se    jack1 jack1.se    jack2     boot  boot.se n
# 1  2014    1      25 41.80000 12.405160 37.00000 6.928203 44.02381 30.33528 4.170934 7
# 2  2014    2      18 20.20408  2.477242 23.14286 2.474358 23.26190 20.78329 1.731019 7
# 3  2014    3      22 43.00000 16.039015 34.00000 8.071113 41.61905 27.18339 4.399938 7
# 4  2014    4      18 39.42857 20.849314 26.57143 5.897076 32.52381 21.67379 3.058999 7
# 5  2014    5      21 26.78571  5.122728 28.71429 4.873921 31.57143 24.71333 2.633834 7
# 6  2014    6      22 34.34286  9.674406 32.28571 6.821335 37.88095 26.63556 3.967503 7
# 7  2014    7      12 13.71429  2.321154 15.42857 2.157096 15.90476 13.78172 1.331006 7
# 8  2014    8      15 18.00000  3.070598 21.00000 3.464102 21.83333 18.06336 2.510976 7
# 9  2014    9      23 30.40816  5.933378 32.42857 5.022399 36.11905 27.46813 2.922924 7
# 10 2014   10      20 41.00000 16.039015 32.00000 5.707138 39.61905 25.16084 2.862035 7
# 11 2014   11      19 43.14286 19.944481 30.14286 6.356742 37.64286 23.74595 3.203750 7
# 12 2014   12      20 34.28571 12.804123 28.57143 4.199125 33.92857 23.78604 2.263036 7
# 13 2014   13      24 48.77143 17.084443 38.57143 7.553888 47.73810 30.27546 3.427311 7
# 14 2014   14      22 34.96429 10.774962 31.42857 4.285714 36.90476 26.20359 2.352035 7
# 15 2014   15      26 67.28571 31.788630 40.57143 9.646825 50.92857 32.12833 4.953348 7
# 16 2014   16      21 28.14286  6.004533 29.57143 3.768830 33.14286 25.01111 2.215218 7
# 17 2014   17      18 27.14286  8.877090 24.85714 3.319700 28.78571 21.10871 1.827654 7
# 18 2014   18      19 44.92857 24.628101 28.42857 6.767268 35.09524 23.03094 3.667259 7
# 19 2014   19      23 37.48571 11.000809 34.14286 5.938460 40.45238 27.94100 3.093979 7
# 20 2014   20      25 33.04762  5.930521 36.14286 5.488392 40.07143 30.31514 3.196420 7
# 21 2014   21      27 37.50000  7.468685 39.00000 6.279217 44.23810 32.59998 3.541380 7
# 22 2014   22      21 49.00000 22.656861 33.00000 6.279217 41.21429 26.10577 2.939562 7
# 23 2014   23      32 49.35714 11.043773 47.42857 7.016017 55.52381 38.96231 3.599812 7
# 24 2014   24      24 56.14286 25.535090 36.85714 6.337449 45.78571 29.48547 3.309005 7




#
########################## 3. Each plot - 2015 ####

raw_int_2015$Round = factor(raw_int_2015$Round)

plot_2015_ls = split(raw_int_2015, raw_int_2015$Plot)

plot_mat_2015_ls = lapply(plot_2015_ls, function(x) {
  mat1 = tapply(x$Record_ID, list(x$Round, x$Specimen_ID), length)
  mat1[is.na(mat1)] = 0
  mat1
})

head(plot_mat_2015_ls[1])
names(plot_mat_2015_ls)
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24"

plot(specaccum(plot_mat_2015_ls[[1]], method="random"))
# cycling through all plots shows none reach asymptote.


# calculate extrapolated richness values:
plot_2015_ex = Reduce(rbind,
                      Map(function(x, y) {
                        df1 = specpool(x)
                        df1$Plot = paste0(y)
                        df1$Year = "2015"
                        rownames(df1) = NULL
                        df1[, c(11, 10, 1:9)]
                      }, plot_mat_2015_ls, names(plot_mat_2015_ls)))
plot_2015_ex
#    Year Plot Species      chao   chao.se    jack1 jack1.se    jack2     boot  boot.se n
# 1  2015    1      18  24.12245  5.145009 26.57143 6.712766 29.54762 22.08311 3.996655 7
# 2  2015    2      17  25.67857  7.809127 24.71429 4.314191 28.76190 20.49855 2.402844 7
# 3  2015    3      20  22.62500  2.724670 26.00000 4.720775 26.23810 23.21791 2.893865 7
# 4  2015    4      22  42.57143 17.398100 32.28571 7.068181 39.07143 26.44864 3.657134 7
# 5  2015    5      17  21.20000  4.219953 23.00000 2.927700 25.02381 19.93597 2.049986 7
# 6  2015    6      18  39.42857 20.849314 26.57143 5.283783 32.52381 21.68852 2.730762 7
# 7  2015    7      15  24.14286  8.877090 21.85714 3.568570 25.78571 18.06361 2.012068 7
# 8  2015    8      20  28.57143  7.256243 28.57143 3.989783 32.73810 23.91374 2.245156 7
# 9  2015    9      17  20.00000  3.070598 23.00000 3.207135 23.83333 20.06602 2.399000 7
# 10 2015   10      17  78.71429 71.693224 27.28571 6.161102 35.26190 21.19657 2.932142 7
# 11 2015   11      20  28.57143  7.256243 28.57143 4.398516 32.73810 23.93583 2.521692 7
# 12 2015   12      12  12.64286  1.123724 14.57143 1.979487 13.14286 13.63138 1.668553 7
# 13 2015   13      21  72.85714 60.982768 30.42857 4.668610 37.69048 24.92415 2.332230 7
# 14 2015   14      18  35.28571 15.017903 27.42857 4.668610 33.50000 22.06628 2.417669 7
# 15 2015   15      15  18.85714  4.182690 20.14286 2.799417 22.04762 17.48145 1.821975 7
# 16 2015   16      16  27.57143 10.757057 23.71429 5.834451 28.35714 19.40368 3.202663 7
# 17 2015   17      18  29.57143 10.757057 25.71429 4.873921 30.35714 21.40899 2.700498 7
# 18 2015   18      12  14.28571  3.041661 15.42857 2.157096 16.50000 13.70660 1.396181 7
# 19 2015   19      18 102.00000 95.686990 30.00000 6.803361 39.40476 22.87890 3.344755 7
# 20 2015   20      15  25.50000 11.341296 21.00000 2.618615 24.80952 17.63180 1.405056 7
# 21 2015   21      16  77.71429 71.693224 26.28571 5.421047 34.26190 20.19907 2.494185 7
# 22 2015   22      24  60.57143 28.579047 37.71429 6.259686 47.35714 29.76352 2.886501 7
# 23 2015   23      24  40.80000 12.405160 36.00000 6.803361 43.02381 29.29815 3.644231 7
# 24 2015   24      23  39.80000 12.405160 35.00000 5.707138 42.02381 28.25336 2.977502 7




####

