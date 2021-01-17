###############################################################################
### R script for running PARTYsocv2.0 and PARTYsocv2.0eu statistical models
###############################################################################
### Written by Lauric Cécillon (lauric.cecillon@inrae.fr)
###############################################################################
### This script accompanies the draft: 
### "Partitioning soil organic carbon into its centennially active and stable fractions 
### with statistical models based on Rock-Eval® thermal analysis (PARTYSOCv2.0 and PARTYSOCv2.0EU)"
### Submitted to Geoscientific Model Development
### By Lauric Cécillon et al.
#####################################################################################################

#####################################################################################################
### STEP 1: load csv data file, create R variables, test correlations, select predictor variables
#####################################################################################################
### Specify pathway for data
setwd("set_your_working_directory")
getwd()
### Empty global environment
rm(list=ls())

### Read csv file
ltbf_re6<-read.csv("PARTYsocv2.0_GMD.csv", header=T, sep=";", row.names=1)
dim(ltbf_re6)
colnames(ltbf_re6)
View(ltbf_re6)

### Potential predictors variables: 40 Rock-Eval parameters
ltbf_re6_predictors<-ltbf_re6[,c(5:44)]
dim(ltbf_re6_predictors)
colnames(ltbf_re6_predictors)
View(ltbf_re6_predictors)

### Response variables = stable SOC proportion and its standard deviation
ltbf_response_variables<-ltbf_re6[,c(48,49)]
dim(ltbf_response_variables)
colnames(ltbf_response_variables)
View(ltbf_response_variables)
rownames(ltbf_response_variables)
summary(ltbf_response_variables)

############################################################################################
### Spearman's correlations between potential predictor variables and stable SOC proportion
### Data not scaled = same spearman rho coefficient than scaled variables
### Data selection (n = 105; with 15 samples per site)
corrrrelations<-matrix(0, nrow = 40, ncol = 2)
rownames(corrrrelations)<-colnames(ltbf_re6_predictors)
colnames(corrrrelations)<-c("rho","p-value")
View(corrrrelations)

for (i in 1:40) {
  corrrrelations[i,1]<-cor.test(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),i],
                                ltbf_response_variables$cpsoc_sample_proportion[which(ltbf_re6$dataset105_7sites_15samples=="yes")],
                                method = c("spearman"))$estimate
  corrrrelations[i,2]<-cor.test(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),i],
                                ltbf_response_variables$cpsoc_sample_proportion[which(ltbf_re6$dataset105_7sites_15samples=="yes")], 
                                method = c("spearman"))$p.value
}

View(corrrrelations)

#####################################################################################################
### Selection of Rock-Eval predictor variables (uses the "corrrrelations" variable created above)
### IMPORTANT!!! Remove RE6 variables with absolute spearman rho < 0.5 with the stable SOC proportion 
#   (in the learning set n = 105)
ltbf_re6_predictors<-ltbf_re6_predictors[,which(abs(corrrrelations[,1])>=0.50)]
dim(ltbf_re6_predictors)
colnames(ltbf_re6_predictors)
#   Only 18 Rock-Eval predictor variables will be used to build PARTYsocv2.0 and PARTYsocv2.0eu models
### Go to STEP 2!

#####################################################################################################
### STEP 2: Prepare the calibration and the validation data sets for PARTYsocv2.0 and PARTYsocv2.0eu
#####################################################################################################

####################################################################################################
### STEP 2 for PARTYsocv2.0 (7 sites, 105 samples), using 18 Rock-Eval predictor variables
####################################################################################################

####################################################################################################
#### Internal validation procedure for PARTYsocv2.0
####################################################################################################
set.seed(1)
### ONLY USE FOR UNCERTAINTY ANALYSIS
x <- c(1:105)
x
randX <- sample(x, length(x), replace=FALSE)
randX

# X calibration set
calre6raw<-ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]
rownames(calre6raw)
dim(calre6raw)

# Y calibration set
(cal_cstable_re6<-ltbf_re6$cpsoc_sample_proportion[which(ltbf_re6$dataset105_7sites_15samples=="yes")])
length(cal_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### Random splitting validation procedure for PARTYsocv2.0
####################################################################################################
#### Repeated 15 times, changing seed from 1 to 15
set.seed(1)
x <- c(1:105)
x
randX <- sample(x, length(x), replace=FALSE)
randX

# X calibration set
calre6raw<-ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),][randX[1:75],]
rownames(calre6raw)
dim(calre6raw)

# X validation set
valre6raw<-ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),][randX[76:105],]
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-ltbf_re6$cpsoc_sample_proportion[which(ltbf_re6$dataset105_7sites_15samples=="yes")][randX[1:75]])
length(cal_cstable_re6)
# Check range of CPsoc values for calibration set
summary(cal_cstable_re6)

# Y validation set
(val_cstable_re6<-ltbf_re6$cpsoc_sample_proportion[which(ltbf_re6$dataset105_7sites_15samples=="yes")][randX[76:105]])
length(val_cstable_re6)
# Check range of CPsoc values for validation set
summary(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
### Leave-one-site-out validation procedure for PARTYsocv2.0
####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for ULTUNA out only (PARTYsocv2.0)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")

rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                  ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                  ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                  ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                  ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                  ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Rothamsted out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Grignon out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Versailles out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")
rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for La Cabana out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Bad Lauchstadt out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Askov out only (PARTYsocv2.0)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="cabana"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
### STEP 2 for PARTYsocv2.0eu (6 sites, 90 samples), using 18 Rock-Eval predictor variables
####################################################################################################

####################################################################################################
#### Internal validation procedure for PARTYsocv2.0eu
####################################################################################################
set.seed(1)
# ONLY USE FOR UNCERTAINTY ANALYSIS
x <- c(1:90)
x
randX <- sample(x, length(x), replace=FALSE)
randX

# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")

rownames(calre6raw)
dim(calre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"))
length(cal_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### Random splitting validation procedure for PARTYsocv2.0eu
####################################################################################################
#### Repeated 15 times, changing seed from 1 to 15
set.seed(1)
x <- c(1:90)
x
randX <- sample(x, length(x), replace=FALSE)
randX

# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")[randX[1:60],]
rownames(calre6raw)
dim(calre6raw)

# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")[randX[61:90],]
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")[randX[1:60]])
length(cal_cstable_re6)
# Check range of CPsoc values for calibration set
summary(cal_cstable_re6)

# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")[randX[61:90]])
length(val_cstable_re6)
# Check range of CPsoc values for validation set
summary(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
### Leave-one-site-out validation procedure for PARTYsocv2.0eu
####################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for ULTUNA out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna")) 
length(val_cstable_re6)
### Go to STEP 3!

######################################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Versailles out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")) 
length(val_cstable_re6)
### Go to STEP 3!

######################################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Askov out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov")) 
length(val_cstable_re6)
### Go to STEP 3!

######################################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Grignon out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")) 
length(val_cstable_re6)
### Go to STEP 3!

######################################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Rothamsted out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted")) 
length(val_cstable_re6)
### Go to STEP 3!

######################################################################################################################
#### ATTENTION!!! Leave-one-site-out: code for Bad Lauchstadt out only (PARTYsocv2.0eu)
# X calibration set
set.seed(1)
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt")) 
length(val_cstable_re6)
### Go to STEP 3!

########################################################################################################
#### Sensitivity of model performance to the reference sites included in the learning set of the model 
########################################################################################################
#### ATTENTION!!! code for Versailles out 
####              with a learning set composed of 4 sites (Grignon, Ultuna, Rothamsted, Askov)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")) 
length(val_cstable_re6)
### Go to STEP 3!

#####################################################################################################################
#### ATTENTION!!! code for Versailles out 
####              with a learning set composed of 3 sites (Grignon, Ultuna, Askov)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! code for Versailles out 
####              with a learning set composed of 2 sites (Grignon, Ultuna)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! code for Grignon out 
####              with a learning set composed of 4 sites (Versailles, Ultuna, Askov, Rothamsted)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! code for Grignon out 
####              with a learning set composed of 3 sites (Versailles, Ultuna, Askov)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")) 
length(val_cstable_re6)
### Go to STEP 3!

####################################################################################################
#### ATTENTION!!! code for Grignon out 
####              with a learning set composed of 2 sites (Versailles, Askov)
set.seed(1)
# X calibration set
calre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),],
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles")


rownames(calre6raw)
dim(calre6raw)
# X validation set
valre6raw<-subset(ltbf_re6_predictors[which(ltbf_re6$dataset105_7sites_15samples=="yes"),], 
                  subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")
dim(valre6raw)
rownames(valre6raw)

# Y calibration set
(cal_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion,
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                           ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"))
length(cal_cstable_re6)
# Y validation set
(val_cstable_re6<-subset(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion, 
                         subset = ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon")) 
length(val_cstable_re6)
### Go to STEP 3!

#####################################################################################################
### STEP 3: Scaling the data sets and running PARTYsocv2.0 and PARTYsocv2.0eu
#####################################################################################################

#####################################################################################################
### STEP 3a: Scaling data
#####################################################################################################
# scaling the calibration set 
ltbf_cal_scaled<-scale(calre6raw, center = TRUE, scale = TRUE)
str(ltbf_cal_scaled)
attributes(ltbf_cal_scaled)$"scaled:center"
attributes(ltbf_cal_scaled)$"scaled:scale"
View(ltbf_cal_scaled)
# new variable without attributes with the scaled calibration dataset
ltbf_calib_scaled<-matrix(ncol = ncol(calre6raw), nrow = nrow(calre6raw), 0)
rownames(ltbf_calib_scaled)<-rownames(calre6raw)
colnames(ltbf_calib_scaled)<-colnames(calre6raw)
for (i in 1:ncol(calre6raw))
{
  ltbf_calib_scaled[,i]<-ltbf_cal_scaled[,i]
}
str(ltbf_calib_scaled)
# how it works
calre6raw[1,1]
(calre6raw[1,1]-attributes(ltbf_cal_scaled)$"scaled:center"[1])/attributes(ltbf_cal_scaled)$"scaled:scale"[1]
ltbf_cal_scaled[1,1]
### (USE ONLY FOR VALIDATION)
# scaling the validation set with the attributes of the calibration dataset
ltbf_val_scaled<-matrix(ncol = ncol(valre6raw), nrow = nrow(valre6raw), 0)
rownames(ltbf_val_scaled)<-rownames(valre6raw)
colnames(ltbf_val_scaled)<-colnames(valre6raw)
for (i in 1:ncol(valre6raw))
{
  ltbf_val_scaled[,i]<-(valre6raw[,i]-attributes(ltbf_cal_scaled)$"scaled:center"[i])/attributes(ltbf_cal_scaled)$"scaled:scale"[i]
}
View(ltbf_val_scaled)
### Go to STEP 3b!

###########################################################################################
### STEP 3b: Running PARTYsocv2.0 and PARTYsocv2.0eu (random forests regression models)
###          (using datasets created at STEP 2 and scaled at STEP 3a)
###########################################################################################
library(randomForest)
### CALIBRATION DATASET (Y et X variables)
# convert data object to data.frame
cal_ltbfRE6_data <- as.data.frame(ltbf_calib_scaled)
dim(cal_ltbfRE6_data)
# final calibration dataset
cal_data<-cbind(cal_cstable_re6,cal_ltbfRE6_data)
dim(cal_data)
colnames(cal_data)

### VALIDATION DATASET (Y et X variables) (ONLY FOR VALIDATION)
# convert data object to data.frame (ONLY FOR VALIDATION)
val_ltbfRE6_data <- as.data.frame(ltbf_val_scaled)
dim(val_ltbfRE6_data)
# final validation dataset (ONLY FOR VALIDATION)
val_data<-cbind(val_cstable_re6,val_ltbfRE6_data)
dim(val_data)
### Go to STEP 3b1!

###########################################################################################
### Step 3b1: random forest regression model on calibration data
rf_model <- randomForest(cal_cstable_re6 ~ ., data = cal_data, ntree=1000, importance=TRUE)
### Variable importance
rf_model$importance
rf_model$importanceSD
varImpPlot(rf_model, scale=FALSE)
### SCALE SHOULD BE SET TO FALSE 
importance(rf_model, type=1, scale=FALSE)
importance(rf_model, type=2, scale=FALSE)
rf_variable_importance<-importance(rf_model,scale=FALSE)
colnames(rf_variable_importance)<-c("%IncMSE","IncNodePurity")
write.csv2(rf_variable_importance,"varimp_partysoc.csv")
### Go to STEP 3b2 (for the UNCERTAINTY ANALYSIS) or directly to STEP 3b5!

###########################################################################################
### Step 3b2: USE ONLY FOR UNCERTAINTY ANALYSIS 
###           (following a methodology adapted from Coulston et al. 2016, see Cécillon et al., 2018, Biogeosciences)
###           bootstrap resampling to parameterize a large number of random forest models
###           Draw nboot bootstraps of the calibration dataset (sample with replacement)
# Matrix with sample indices selected in each bootstrap
(b<-c(1:nrow(cal_data)))
nboot<-2000
(bootstraps<-matrix(0,nboot,nrow(cal_data)))
for (i in 1:nboot) 
{
  bootstraps[i,]<-sample(b,length(b), replace=TRUE)
}
dim(bootstraps)

# Find the frequency of each indice in each bootstrap (and non selected indices within each bootstrap) 
library(boot)
(mat_freq<-freq.array(bootstraps))

# Store "out-of-bag" indices of each bootstrap in a list
out_of_bag<-list()
for (i in 1:nboot)
{
  out_of_bag[[i]]<-which(mat_freq[i,]==0)
}

########
# Create a list of the bootstrapped calibration data (Y=cstable) & (X=RE6 variables) 
# In each bootstrap of Y, we put a random value of Y generated by the rnorm function 
# with a mean calculated value and the estimated standard deviation

##############################
# ATTENTION!!! Use only for model PARTYsocv2.0 (7 sites)
Y<-matrix(nrow=105, ncol=nboot, 0)
rownames(Y)<-rownames(ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),])
for (i in 1:105)
{
  Y[i,]<-rnorm(nboot, mean = ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$cpsoc_sample_proportion[i], sd = ltbf_response_variables[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$sd_cpsoc_sample_proportion[i])
}
View(Y)
Y[36,]
randX[1:105]
### Y_cal for model PARTYsocv2.0
Y_cal<-Y[randX[1:105],]

##############################
### ATTENTION!!! Use only for model PARTYsocv2.0eu (6 sites)
Y<-matrix(nrow=90, ncol=nboot, 0)
rownames(Y)<-rownames(ltbf_response_variables[which(ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"),])
for (i in 1:90)
{
  Y[i,]<-rnorm(nboot, 
               mean = ltbf_response_variables[which(ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                                                      ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"),]$cpsoc_sample_proportion[i], 
               sd = ltbf_response_variables[which(ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="versailles"|
                                                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="ultuna"|
                                                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="askov"|
                                                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="rothamsted"|
                                                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="bad_lauchstadt"|
                                                    ltbf_re6[which(ltbf_re6$dataset105_7sites_15samples=="yes"),]$site=="grignon"),]$sd_cpsoc_sample_proportion[i])
}
View(Y)
Y[36,]

### ATTENTION!!! Y_cal for PARTYsocv2.0eu
Y_cal<-Y[randX[1:90],]
###############################

boot_calibration<-list()
for (i in 1:nboot)
{
  boot_calibration[[i]]<-cbind(Y[bootstraps[i,],i],cal_data[bootstraps[i,],2:ncol(cal_data)])
  # Make replicated sample names unique for future use in the random forests (RF) model
  rownames(boot_calibration[[i]])<-make.unique(rownames(boot_calibration[[i]]), sep = "_")
  # RF model needs to find Y variable name, set this name, same for all bootstraps 
  colnames(boot_calibration[[i]])<-c("y_cal_boot",colnames(cal_data[,2:ncol(cal_data)]))
}

dim(boot_calibration[[1]])
boot_calibration[[1]]
rownames(boot_calibration[[1]])
colnames(boot_calibration[[1]])

# Build a list of Rock-Eval data for validation samples (X variables only) 
# For each bootstrap (out of bag), the size of out of bag vectors can cary (ca. 37%)
out_boot_prediction<-list()
for (i in 1:nboot)
{
  out_boot_prediction[[i]]<-cal_data[out_of_bag[[i]],2:ncol(cal_data)]
}

dim(out_boot_prediction[[1]])
dim(out_boot_prediction[[2]])

# Build a list of observed stable SOC fraction data in validation samples (Y variable only)
# for each bootstrap (out of bag) 
# In each bootstrap, put a random value for Y generated by the rnorm function
# with mean the calculated value and standard deviation the estimaded error
out_boot_actual_c_ctable<-list()
for (i in 1:nboot)
{
  out_boot_actual_c_ctable[[i]]<-Y[out_of_bag[[i]],i]
}

length(out_boot_actual_c_ctable[[1]])
length(out_boot_actual_c_ctable[[2]])

### Make nboot random forest models (of 1000 trees) for each of the nboot boostraps
# In each bootstrap, we put a random value for Y
# We store the results of the nboot RF models in a list
?randomForest
rf<-list()
for (i in 1:nboot)
{
  rf[[i]] <- randomForest(y_cal_boot ~ ., data=boot_calibration[[i]],
                          ntree=1000,
                          importance=TRUE)
}

names(rf[[1]])
mean(rf[[1]]$rsq)
mean(rf[[1]]$mse)
varImpPlot(rf[[1]], scale=FALSE)
varImpPlot(rf[[2]], scale=FALSE)
varImpPlot(rf[[3]], scale=FALSE)

### We make a prediction of the stable SOC proportion for the nboot set of out of bag samples
# Store the predictions of the nboot RF models in a list
rf_pred<-list()
for (i in 1:nboot)
{
  #with predict.all = TRUE, keep the prediction of each tree for each prediction 
  rf_pred[[i]] <- predict(rf[[i]], out_boot_prediction[[i]], predict.all=TRUE)
}

names(rf_pred[[1]])
rf_pred[[1]]$individual[1,]
#mean of prédictions of the 1000 trees (individuals) for the first out-of-bag sample
mean(rf_pred[[1]]$individual[1,])
#same result of the aggregate value
rf_pred[[1]]$aggregate[1]
dim(rf_pred[[1]]$individual)

# For each out-of-bag sample (varying number) of each bootstrap (nboot), 
# Store the value of the variance of predictions of the 1000 trees in a list
# Use the function apply to apply the function "var" on the rowas of the data frame (parameter=1)
variance<-list()
for (i in 1:nboot)
{
  variance[[i]]<-apply(rf_pred[[i]]$individual,1,var)
}

variance[[1]]
### Go to STEP 3b3!

################################################################################
### Step 3b3: USE ONLY FOR UNCERTAINTY ANALYSIS 
###           (following a methodology adapted from Coulston et al. 2016, see Cécillon et al., 2018, Biogeosciences)
# Store for each out-of-bag sample of each bootstrap:
# - observed values of stable SOC proportion (out_boot_actual_c_ctable), 
# - predicted values (rf_pred$aggregate) 
# - variance of predicted values on the 1000 trees 
# "Step 3b3 yields an error assessment dataset"
error_assessment_dataset<-list()
for (i in 1:nboot)
{
  error_assessment_dataset[[i]]<-cbind(out_boot_actual_c_ctable[[i]], rf_pred[[i]]$aggregate, variance[[i]])
}

error_assessment_dataset[[1]]
error_assessment_dataset[[2]]
### Go to STEP 3b4!

################################################################################
### Step 3b4: USE ONLY FOR UNCERTAINTY ANALYSIS 
###           (following a methodology adapted from Coulston et al. 2016, see Cécillon et al., 2018, Biogeosciences)
# Compute and store a t value for each out-of-bag sample of each bootstrap.
# t is calculated according to the formula t = sqrt((yobs-ypred)²/var(y))
to<-list()
for (i in 1:nboot)
{
  to[[i]]<-sqrt(((error_assessment_dataset[[i]][,1]-error_assessment_dataset[[i]][,2])^2)/error_assessment_dataset[[i]][,3])
}

to[[1]]
to[[2]]

### Cumulative distribution of t, and estimation of t (95%) with the MONTE CARLO approach
to_cumul_MC<-vector()
for (i in 1:nboot)
{
  to_cumul_MC<-c(to_cumul_MC,to[[i]])
}

ecdf(to_cumul_MC)
plot(ecdf(to_cumul_MC), verticals = TRUE, las=1, lwd=6, xlim=c(0,6))
segments(quantile(to_cumul_MC, probs = 0.95), 0.95, quantile(to_cumul_MC, probs = 0.95), -0.5, col = "gray60", lty=3, lwd=6)
segments(quantile(to_cumul_MC, probs = 0.95), 0.95, -1.7, 0.95, col = "gray60", lty=3, lwd=6)

summary(to_cumul_MC)
# here is the value of t (95%) with the MONTE CARLO approach
(t_value_MC<-quantile(to_cumul_MC, probs = 0.95))

### Cumulative distribution of t, and estimation of t (95%) with the BOOTSTRAP approach
to_cumul_boot<-vector()
for (i in 1:nboot)
{
  to_cumul_boot[i]<-quantile(to[[i]], probs = 0.95)
}

# here is the value of t (95%) with the BOOTSTRAP method
(t_value_BOOT<-mean(to_cumul_boot))
### Go to STEP 3b5!

#####################################################################################################
### Step 3b5: (USE ONLY FOR VALIDATION) 
###           make a prediction of the stable SOC proportion, including prediction error, 
###           for a new observation (on the validation set)
rf_model_pred <- predict(rf_model, val_data[,2:ncol(val_data)], predict.all=TRUE)

names(rf_model_pred)
dim(rf_model_pred$individual)
# store standard deviations of predictions for each validation sample (over 1000 trees)
sd_pred_model<-vector()
for (i in 1:nrow(rf_model_pred$individual))
{
  sd_pred_model[i]<-sd(rf_model_pred$individual[i,])
}

sd_pred_model
### Go to STEP 4!

###################################################################################
### STEP 4: Compute model performance statistics (see definitions in the paper)
###################################################################################
### R²oob
mean(rf_model$rsq)
### RMSEPoob
(RMSEC<-sqrt((sum(di2<-(cal_data[,1]-rf_model$predicted)^2))/length(cal_data[,1])))
### R²
cor(rf_model_pred$aggregate,val_data[,1], method = c("pearson"))^2
### RMSEP
(RMSEP<-sqrt((sum(di2<-(val_data[,1]-rf_model_pred$aggregate)^2))/length(val_data[,1])))
### rRMSEP
(RMSEP/mean(val_data[,1]))
### RPIQ 
(RPIQ<-(quantile(val_data[,1],0.75)-quantile(val_data[,1],0.25))/RMSEP)
### BIAS
(BIAS<-mean(rf_model_pred$aggregate)-mean(val_data[,1]))

########################################################################################
### Get results of the internal validation procedure
(validation<-cbind(rf_model$predicted,cal_data[,1]))
colnames(validation)<-c("prediction", "reference_value")
View(validation)

########################################################################################
### Get results of the leave-one-site-out or of the random sampling validation procedure
(validation<-cbind(rf_model_pred$aggregate,val_data[,1]))
colnames(validation)<-c("prediction", "reference_value")
View(validation)
