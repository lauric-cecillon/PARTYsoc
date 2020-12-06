################################################################################
#   PARTYsocv2.0eu
#   A statistical model partitioning soil organic carbon into its
#   centennially active and stable fractions based on Rock-Eval thermal analysis 
#   R script by Lauric Cécillon 
#   Contact: lauric.cecillon@inrae.fr
################################################################################

### Specify pathway for data
setwd("HERE specify your working directory")
rm(list=ls())

### ATTENTION!!! Use only on csv file with 18 Rock-Eval parameters
### ATTENTION: S2 should be in gC_kg 
#   (from home made R script, not from the Vinci Technologies software where it is in gCH/kg)
### PC (mgC/g), HI, OI RE6, TOC_RE6 (gC/gkg) are derived from Vinci Technologies software
#   ordered like this:

# "t_ch_pyr_70" "t_ch_pyr_90"    
# "t_co2_pyr_30" "t_co2_pyr_50" "t_co2_pyr_70" "t_co2_pyr_90"
# "t_co_ox_70" 
# "t_co2_ox_50" "t_co2_ox_70" "t_co2_ox_90"
# "pseudos1_gc_kg" "s2_gc_kg" "s2_pc"
# "hi" "hi_oi_re6" "pc_g_kg" "pc_toc_re6" "toc_re6_g_kg"  

re6_predictors<-read.csv("name_of_file_with_18_Rock-Eval_predictor_variables.csv", 
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
  new_data_scaled[,i]<-(re6_predictors[,i]-attributes(ltbf_cal_scaled_6sitesEU)$"scaled:center"[i])/
                        attributes(ltbf_cal_scaled_6sitesEU)$"scaled:scale"[i]
}
View(new_data_scaled)
dim(new_data_scaled)

###########################################################################################
### Random Forests regression model
### Make a prediction of the centennially stable SOC proportion, including error, 
#  for a new soil sample
### NEW DATASET FOR RF model (X variables)
library(randomForest)
rf_model_6sitesEU_pred <- predict(rf_model_6sitesEU, new_data_scaled, predict.all=TRUE)
names(rf_model_6sitesEU_pred)
rf_model_6sitesEU_pred$aggregate
dim(rf_model_6sitesEU_pred$individual)

# store standard deviations of predictions for each new sample (over 1000 trees)
sd_pred_model<-vector()

for (i in 1:nrow(rf_model_6sitesEU_pred$individual))
{
  sd_pred_model[i]<-sd(rf_model_6sitesEU_pred$individual[i,])
}

sd_pred_model

# store prediction and compute final error for each validation sample
final_pred_error<-matrix(nrow=nrow(rf_model_6sitesEU_pred$individual), ncol=3,0)
# Modèle PARTYsocv2.0eu
colnames(final_pred_error)<-c("CPSOC_proportion_PARTYsocv2.0eu", "CPSOC_proportion_PARTYsocv2.0eu_ci95_mc", "CPSOC_proportion_PARTYsocv2.0eu_ci95_boot")

rownames(final_pred_error)<-names(rf_model_6sitesEU_pred$aggregate)
final_pred_error[,1]<-rf_model_6sitesEU_pred$aggregate
final_pred_error[,2]<-sd_pred_model*t_value_MC_6sitesEU
final_pred_error[,3]<-sd_pred_model*t_value_BOOT_6sitesEU
View(final_pred_error)

write.csv2(final_pred_error, "name_of_result_file.csv")
