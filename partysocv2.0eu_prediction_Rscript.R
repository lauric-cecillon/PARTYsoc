###############################################################################
### R script by Lauric Cécillon (lauric.cecillon@inrae.fr)
###############################################################################

### Specify pathway for data
### Specify pathway for data
setwd("HERE specify your working directory")
getwd()
rm(list=ls())

### ATTENTION!!! Use only on csv file with 18 parameters 
### ATTENTION: S2 should be in gC_kg (from our script, not from Vinci software where it is in gCH/kg!)
### PC (mgC/g), HI, OI RE6, TOC_RE6 (gC/gkg) are derived from Vinci software
#   ordered like this:

# "t_ch_pyr_70" "t_ch_pyr_90"    
# "t_co2_pyr_30" "t_co2_pyr_50" "t_co2_pyr_70" "t_co2_pyr_90"
# "t_co_ox_70" 
# "t_co2_ox_50" "t_co2_ox_70" "t_co2_ox_90"
# "pseudos1_gc_kg" "s2_gc_kg" "s2_pc"
# "hi" "hi_oi_re6" "pc_g_kg" "pc_toc_re6" "toc_re6_g_kg"  

re6_predictors<-read.csv("name_of_file_with_18_Rock-Eval_predictor_variables.csv", header=T, sep=";", row.names=1)
dim(re6_predictors)
colnames(re6_predictors)
rownames(re6_predictors)
View(re6_predictors)

# OPEN R.data file "PARTYsocv2.0_PARTYsocv2.0eu_key_data_for_prediction"

#############################################################################
### Scaling data
# scaling the new sample set with the attributes of the calibration dataset
new_data_scaled<-matrix(ncol = ncol(re6_predictors), nrow = nrow(re6_predictors), 0)
rownames(new_data_scaled)<-rownames(re6_predictors)
colnames(new_data_scaled)<-colnames(re6_predictors)
for (i in 1:ncol(re6_predictors))
{
  new_data_scaled[,i]<-(re6_predictors[,i]-attributes(ltbf_cal_scaled_6sitesEU)$"scaled:center"[i])/attributes(ltbf_cal_scaled_6sitesEU)$"scaled:scale"[i]
}
View(new_data_scaled)
dim(new_data_scaled)

#########################################
# Mahalanobis distance

########## MODELE 6 SITES GRI VERS ULT ASK ROTH BAD
# Travail sur sélection de variables les plus importantes dans le modèle
# PARTYsocv2.0eu
# les 5 variables les plus importantes sont :
# S2, T50CO2PYR, PC, S2/PC, HI/OIRE6
colnames(ltbf_cal_scaled_6sitesEU)
(moyx<-colMeans(ltbf_cal_scaled_6sitesEU[,c(12,4,16,13,15)]))
(Sx<-cov(ltbf_cal_scaled_6sitesEU[,c(12,4,16,13,15)]))
(D2<-mahalanobis(ltbf_cal_scaled_6sitesEU[,c(12,4,16,13,15)],moyx,Sx))
hist(D2)
mean(D2)
# Distance max au sein du jeu de données LTBF
max(D2)

# Distance des nouvelles données par rapport au jeu de données LTBF 
# (basée sur les 5 variables les plus importantes du modèle)
(Distance_new_data<-mahalanobis(new_data_scaled[,c(12,4,16,13,15)],moyx,Sx))
hist(Distance_new_data)
mean(Distance_new_data)
max(Distance_new_data)
write.csv2(Distance_new_data, "maha_distances_songchao_ltbf_aial_proterr_max_14_45234.csv")

# Sélection des échantillons sur le critère 
# 2 fois la distance max au sein du jeu de données LTBF
SELECTED_new_data_scaled<-subset(new_data_scaled, subset = Distance_new_data<=2*max(D2))
dim(SELECTED_new_data_scaled)
rownames(SELECTED_new_data_scaled)


###########################################################################################
### RF MODEL
### Step 5: make a prediction, including error, for a new observation (on the new data set)
### NEW DATASET FOR RF model (X variables)
# convert data object to data.frame
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

