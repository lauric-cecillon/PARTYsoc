#####################################################################################################
### R script for running PARTYsocv2.0eu statistical model
### "A statistical model partitioning soil organic carbon into its
### centennially active and stable fractions, based on Rock-Eval thermal analysis"
#####################################################################################################
### Written by Lauric Cécillon (lauric.cecillon@inrae.fr)
#####################################################################################################
### This script accompanies the draft: 
### "Partitioning soil organic carbon into its centennially active and stable fractions 
### with statistical models based on Rock-Eval® thermal analysis (PARTYSOCv2.0 and PARTYSOCv2.0EU)"
### Submitted to Geoscientific Model Development
### By Lauric Cécillon et al.
#####################################################################################################

### Specify pathway for data
setwd("set_your_working_directory")
rm(list=ls())

### ATTENTION!!! Use only on csv file with 18 Rock-Eval parameters 
### calculated as described in the above-mentionned draft
### Parameters should be named as described below, and ordered accordingly, one variable per column
# "T70_HC_PYR"	"T90_HC_PYR"	
# "T30_CO2_PYR"	"T50_CO2_PYR"	"T70_CO2_PYR"	"T90_CO2_PYR"
# "T70_CO_OX"	
# "T50_CO2_OX"	"T70_CO2_OX"	"T90_CO2_OX"	
# "PseudoS1"	"S2"	"S2_PC"	"HI"	"HI_OIre6"	"PC"	"PC_TOCre6"	"TOCre6"

# Load csv file with 18 Rock-Eval parameters (predictor variables) of new topsoil samples
re6_predictors<-read.csv("filename_18_rockeval_predictor_variables.csv", 
                         header=T, sep=";", row.names=1)
dim(re6_predictors)
colnames(re6_predictors)
rownames(re6_predictors)
View(re6_predictors)

### Load the Rdata file provided that contains all necessary information 
#   to run the PARTYsocv2.0eu statistical model
load("PARTYsocv2.0eu_key_data_for_prediction.Rdata")

#############################################################################
### Scaling the 18 Rock-Eval predictor variables of the new soil sample set
#   using the scaling attributes of the PARTYsocv2.0eu learning set (n = 90)
new_data_scaled<-matrix(ncol = ncol(re6_predictors), nrow = nrow(re6_predictors), 0)
rownames(new_data_scaled)<-rownames(re6_predictors)
colnames(new_data_scaled)<-colnames(re6_predictors)
for (i in 1:ncol(re6_predictors))
{
  new_data_scaled[,i]<-(re6_predictors[,i]-attributes(ltbf_cal_scaled)$"scaled:center"[i])/
    attributes(ltbf_cal_scaled)$"scaled:scale"[i]
}
View(new_data_scaled)
dim(new_data_scaled)

###########################################################################################
### Random Forests regression model
### Make a prediction of the centennially stable SOC proportion, including error, 
### for a new soil sample
library(randomForest)
partysocv2.0eu_pred <- predict(rf_model, new_data_scaled, predict.all=TRUE)
names(partysocv2.0eu_pred)
partysocv2.0eu_pred$aggregate
dim(partysocv2.0eu_pred$individual)

# Store standard deviations of predictions for each new sample (over 1000 trees)
sd_pred_model<-vector()
for (i in 1:nrow(partysocv2.0eu_pred$individual))
{
  sd_pred_model[i]<-sd(partysocv2.0eu_pred$individual[i,])
}
sd_pred_model

### Store PARTYsocv2.0eu predictions of the centennially stable SOC proportion in new topsoils 
#   and calculate prediction error (confidence interval; ci95) for each sample, using a Monte-Carlo (mc) or a bootstrap (boot) approach
#   Also calculate the content of the centennially stable SOC fraction, 
#   and the proportion and content of the centennially active SOC fraction
final_pred_error<-matrix(nrow=nrow(partysocv2.0eu_pred$individual), ncol=6,0)
colnames(final_pred_error)<-c("Stable_SOC_proportion_PARTYsocv2.0eu", "Stable_SOC_proportion_PARTYsocv2.0eu_ci95_mc", "Stable_SOC_proportion_PARTYsocv2.0eu_ci95_boot", "Stable_SOC_content_PARTYsocv2.0eu", "Active_SOC_proportion_PARTYsocv2.0eu", "Active_SOC_content_PARTYsocv2.0eu")
rownames(final_pred_error)<-names(partysocv2.0eu_pred$aggregate)
final_pred_error[,1]<-partysocv2.0eu_pred$aggregate
final_pred_error[,2]<-sd_pred_model*t_value_MC
final_pred_error[,3]<-sd_pred_model*t_value_BOOT
final_pred_error[,4]<-partysocv2.0eu_pred$aggregate*re6_predictors$TOCre6
final_pred_error[,5]<-(1-partysocv2.0eu_pred$aggregate)
final_pred_error[,6]<-(1-partysocv2.0eu_pred$aggregate)*re6_predictors$TOCre6
View(final_pred_error)
# Save results as a csv file
write.csv2(final_pred_error, "name_of_result_file.csv")
