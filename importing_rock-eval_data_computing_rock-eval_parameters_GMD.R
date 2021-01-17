#####################################################################################################
### R script for:
### (1) extracting Rock-Eval thermograms and parameters from raw Rock-Eval 6 R.00 and S.00 files 
### (2) computing Rock-Eval parameters (that are not provided by the Vinci Technologies Software)
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


######################################################################################################
### PART 1 : extracting Rock-Eval thermograms and parameters from raw Rock-Eval 6 R.00 and S.00 files
######################################################################################################

############################################################################################################
### Importing baselines of thermograms from raw Rock-Eval 6 R.00 files
############################################################################################################
ls1 <- list.files()
ls1
length(ls1)

library(stringr)
names_samples<-str_replace_all(ls1,".R00","")
names_samples

### Read all R00 files to extract baseline information on pyrolysis and oxidation stage 
### Correction bug identifié pour certains fichiers en mettant l'argument "blank.lines.skip=F"
data_r00 <- lapply(ls1, function(x) read.table(x, header=F, sep="", blank.lines.skip=F))
str(data_r00, 1)

dim(data_r00[[1]])
data_r00[[3]]

# Detect the index of various baselines in each file
mat_pos<-matrix(0,length(names_samples),10)
rownames(mat_pos)<-names_samples
colnames(mat_pos)<-c("ch_pyr_auto","co_pyr_auto","co2_pyr_auto","ch_pyr_manu","co_pyr_manu","co2_pyr_manu","co_ox_auto","co2_ox_auto","co_ox_manu","co2_ox_manu")
View(mat_pos)
for (i in 1:length(names_samples))
{
  mat_pos[i,]<-which(str_detect(data_r00[[i]][,1], "Base"))
}

View(mat_pos)

### Baseline HC_pyr
base_ch_pyr<-matrix(0,length(names_samples),1)
rownames(base_ch_pyr)<-names_samples

for (i in 1:length(names_samples))
{
  base_ch_pyr[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,4]],"Base=",""))
}

base_ch_pyr
write.csv2(base_ch_pyr, "base_ch_pyr.csv")

### Baseline CO_pyr
base_co_pyr<-matrix(0,length(names_samples),1)
rownames(base_co_pyr)<-names_samples

for (i in 1:length(names_samples))
{
  base_co_pyr[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,5]],"Base=",""))
}

base_co_pyr
write.csv2(base_co_pyr, "base_co_pyr.csv")

### Baseline CO2_pyr
base_co2_pyr<-matrix(0,length(names_samples),1)
rownames(base_co2_pyr)<-names_samples

for (i in 1:length(names_samples))
{
  base_co2_pyr[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,6]],"Base=",""))
}

base_co2_pyr
write.csv2(base_co2_pyr, "base_co2_pyr.csv")

### Baseline Co_ox
base_co_ox<-matrix(0,length(names_samples),1)
rownames(base_co_ox)<-names_samples

for (i in 1:length(names_samples))
{
  base_co_ox[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,9]],"Base=",""))
}

base_co_ox
write.csv2(base_co_ox, "base_co_ox.csv")

### Baseline CO2_ox
base_co2_ox<-matrix(0,length(names_samples),1)
rownames(base_co2_ox)<-names_samples

for (i in 1:length(names_samples))
{
  base_co2_ox[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,10]],"Base=",""))
}

base_co2_ox
write.csv2(base_co2_ox, "base_co2_ox.csv")

############################################################################################################
### Importing KFID parameter from raw Rock-Eval 6 R.00 files
############################################################################################################
ls1 <- list.files()
ls1
length(ls1)

library(stringr)
names_samples<-str_replace_all(ls1,".R00","")
names_samples

### Read all R00 files to extract baseline information on pyrolysis and oxidation stage 
### Correction bug identifié pour certains fichiers en mettant l'argument "blank.lines.skip=F"
data_r00 <- lapply(ls1, function(x) read.table(x, header=F, sep="", blank.lines.skip=F))
str(data_r00, 1)

#dim(data_r00[[50]])
#data_r00[[50]]

# Detect the index of kfid in each file
mat_pos<-matrix(0,length(names_samples),1)
rownames(mat_pos)<-names_samples
colnames(mat_pos)<-c("pos_kFid")
View(mat_pos)
for (i in 1:length(names_samples))
{
  mat_pos[i,1]<-which(str_detect(data_r00[[i]][,1], "KFid"))[1]
}

View(mat_pos)

### Extract kFID value for each file
kfid<-matrix(0,length(names_samples),1)
rownames(kfid)<-names_samples

for (i in 1:length(names_samples))
{
  kfid[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,1]],"KFid=",""))
}

kfid
write.csv2(kfid, "kFID.csv")

############################################################################################################
### Importing sample weight from raw Rock-Eval 6 R.00 files
############################################################################################################
ls1 <- list.files()
ls1
length(ls1)

library(stringr)
names_samples<-str_replace_all(ls1,".R00","")
names_samples

### Read all R00 files to extract baseline information on pyrolysis and oxidation stage 
### Correction bug identifié pour certains fichiers en mettant l'argument "blank.lines.skip=F"
data_r00 <- lapply(ls1, function(x) read.table(x, header=F, sep="", blank.lines.skip=F))
str(data_r00, 1)

#dim(data_r00[[50]])
#data_r00[[50]]

# Detect the index of sample weight in each file
mat_pos<-matrix(0,length(names_samples),1)
rownames(mat_pos)<-names_samples
colnames(mat_pos)<-c("pos_sample_weight")
View(mat_pos)
for (i in 1:length(names_samples))
{
  mat_pos[i,1]<-which(str_detect(data_r00[[i]][,1], "Quant"))[1]
}

View(mat_pos)

### Extract sample weight value for each file
sample_weight<-matrix(0,length(names_samples),1)
rownames(sample_weight)<-names_samples

for (i in 1:length(names_samples))
{
  sample_weight[i]<-as.numeric(str_replace_all(data_r00[[i]][,1][mat_pos[i,1]],"Quant=",""))
}

sample_weight
write.csv2(sample_weight, "sample_weight.csv")

############################################################################################################
### Importing Rock-Eval thermograms from raw Rock-Eval 6 S.00 files
############################################################################################################
###########################################
# NB :  First check for potential problems 
#       in the size of the s00 files 
#       (small size may indicate problems)
###########################################
ls1 <- list.files()
ls1
length(ls1)

library(stringr)
names_samples<-str_replace_all(ls1,".S00","")
names_samples

### Read all S00 files to extract information on pyrolysis stage (time, temperature, ch, co, co2)
#   NB: It may bug when S00 file has missing data for several time measurements (a few seconds may be missing)
#   If this is the case, try to reduce the nrows parameter 
#   and identify the sample with bad acquisition using line 25 = time_pyr[nrow(time_pyr),], 
#   the last row of the time table may help in this process! (isolate the sample with a different number)
#   NB2: When only a limited number of seconds are missing, one may reconstruct the signal 
#   by averaging the previous and following temperature and signal intensities for all missing seconds
#   then the modified file should be saved as a text file and stored and treated separately (see line 87)
pyr <- lapply(ls1, function(x) read.table(x, header=F, sep="", skip=2, nrows=1264))
#   For use when problems appear
#   pyr <- lapply(ls1, function(x) read.table(x, header=F, sep="", skip=2, nrows=1261))

tot_pyr <- do.call("cbind", pyr)
View(tot_pyr)
time_pyr<-tot_pyr[,c(T,F,F,F,F)]
View(time_pyr)
#   Check if all last measurement time of the pyrolysis stage is the same for all samples (usually = 1263)
time_pyr[nrow(time_pyr),]
min(time_pyr[nrow(time_pyr),])
max(time_pyr[nrow(time_pyr),])
#   which(time_pyr[nrow(time_pyr),]==1263)
#   names_samples[45]

temp_pyr<-tot_pyr[,c(F,T,F,F,F)]
View(temp_pyr)
### Visual check of the heating ramp (pyrolysis stage) for all samples
#    Identify spectra with problems during the pyrolysis stage
#    (click on the spectra in the plotting window)
#    Returns the indice of the spectra ("ispc")
library(hyperSpec)
tempPYR<-new("hyperSpec", spc = t(temp_pyr))
tempPYR@wavelength
tempPYR@label$.wavelength<-expression(Time /s)
tempPYR@label$spc<-"Temperature /°C"
tempPYR
spc.identify(plot(tempPYR, axis.args = list (las = 1, cex.axis = 1), lines.args = list (pch = 20, lwd = 2), col = "darkblue", spc.nmax = length(names_samples), wl.range = c (min ~ max)), cex = 1.5)
#   Display the name of the sample with problems during the pyrolysis stage 
#   (modify with the corresponding "ispc" number)
names_samples[57]

ch_pyr<-tot_pyr[,c(F,F,T,F,F)]
View(ch_pyr)
co_pyr<-tot_pyr[,c(F,F,F,T,F)]
View(co_pyr)
co2_pyr<-tot_pyr[,c(F,F,F,F,T)]
View(co2_pyr)

### Read all S00 files to extract information on oxidation stage (time, temperature, co, co2)
#   NB: It may bug when S00 file has missing data for several time measurements (a few seconds may be missing)
#   If this is the case, try to reduce the nrows parameter 
#   and identify the sample with bad acquisition using line 57 = time_ox[nrow(time_ox),], 
#   the last row of the time table may help in this process! (isolate the sample with a different number)
#   NB2: When only a limited number of seconds are missing, one may reconstruct the signal 
#   by averaging the previous and following temperature and signal intensities for all missing seconds
#   then the modified file should be saved as a text file and stored and treated separately (see line 87)
ox <- lapply(ls1, function(x) read.table(x, header=F, sep="", skip=1268, nrows=2311))
#   For use when problems appear
#   ox <- lapply(ls1, function(x) read.table(x, header=F, sep="", skip=1268, nrows=2308))

tot_ox <- do.call("cbind", ox)
View(tot_ox)
time_ox<-tot_ox[,c(T,F,F,F)]
View(time_ox)
#   Check if all last measurement time of the oxidation stage is the same for all samples (usually = 2310)
time_ox[nrow(time_ox),]
min(time_ox[nrow(time_ox),])
max(time_ox[nrow(time_ox),])
#   which(time_ox[nrow(time_ox),]==2310)
#   names_samples[45]

temp_ox<-tot_ox[,c(F,T,F,F)]
View(temp_ox)
### Visual check of the heating ramp (oxidation stage) for all samples
#    Identify spectra with problems during the pyrolysis stage
#    (click on the spectra in the plotting window)
#    Returns the indice of the spectra ("ispc")
tempOX<-new("hyperSpec", spc = t(temp_ox))
tempOX@wavelength
tempOX@label$.wavelength<-expression(Time /s)
tempOX@label$spc<-"Temperature /°C"
tempOX
spc.identify(plot(tempOX, axis.args = list (las = 1, cex.axis = 1), lines.args = list (pch = 20, lwd = 2), col = "darkred", spc.nmax = length(names_samples), wl.range = c (min ~ max)), cex = 1.5)
#   Display the name of the sample with problems during the oxidation stage 
#   (modify with the corresponding "ispc" number)
names_samples[74]

co_ox<-tot_ox[,c(F,F,T,F)]
View(co_ox)
co2_ox<-tot_ox[,c(F,F,F,T)]
View(co2_ox)

### Add names_samples as colnames for all data frames
colnames(time_pyr)<-names_samples
colnames(temp_pyr)<-names_samples
colnames(ch_pyr)<-names_samples
colnames(co_pyr)<-names_samples
colnames(co2_pyr)<-names_samples
colnames(time_ox)<-names_samples
colnames(temp_ox)<-names_samples
colnames(co_ox)<-names_samples
colnames(co2_ox)<-names_samples

###########################################################################################################
### Only use for sample with problem during acquisition, after modifying it and saving it as a text file
### Read problematic file (now add the modified text file in the folder)
sample_pyr<-read.table("sample_modif.txt", header=F, sep="", skip=1, nrows=1264)
colnames(sample_pyr)<-c("sample","sample","sample","sample","sample")
View(sample_pyr)
sample_ox<-read.table("sample_modif.txt", header=F, sep="", skip=1267, nrows=2311)
colnames(sample_ox)<-c("sample","sample","sample","sample")
View(sample_ox)

### When necessary, Add the modified file in the data frames
time_pyr<-cbind(time_pyr,sample_pyr[,1])
temp_pyr<-cbind(temp_pyr,sample_pyr[,2])
ch_pyr<-cbind(ch_pyr,sample_pyr[,3])
co_pyr<-cbind(co_pyr,sample_pyr[,4])
co2_pyr<-cbind(co2_pyr,sample_pyr[,5])

time_ox<-cbind(time_ox,sample_ox[,1])
temp_ox<-cbind(temp_ox,sample_ox[,2])
co_ox<-cbind(co_ox,sample_ox[,3])
co2_ox<-cbind(co2_ox,sample_ox[,4])

#######################################################################################################
### Save the 7 CSV files that will be later used to compute RE6 indices
View(time_pyr)
View(temp_pyr)
write.csv2(temp_pyr, "temp_pyr.csv")
View(ch_pyr)
write.csv2(ch_pyr, "ch_pyr.csv")
View(co_pyr)
write.csv2(co_pyr, "co_pyr.csv")
View(co2_pyr)
write.csv2(co2_pyr, "co2_pyr.csv")

View(time_ox)
View(temp_ox)
write.csv2(temp_ox, "temp_ox.csv")
View(co_ox)
write.csv2(co_ox, "co_ox.csv")
View(co2_ox)
write.csv2(co2_ox, "co2_ox.csv")


#########################################################################################
### PART 2 : computing Rock-Eval parameters
#########################################################################################

############################################################################################################
### Compute temperature parameters of the HC_PYR thermogram
############################################################################################################
rm(list=ls())
### Read files
# HC_PYR thermogram
ch_pyr<-read.csv("ch_pyr.csv", header=T, row.names=1, sep=";")
dim(ch_pyr)
View(ch_pyr)
colnames(ch_pyr)
# Oven temperature during pyrolysis
temp_pyr<-read.csv("temp_pyr.csv", header=T, row.names=1, sep=";")
dim(temp_pyr)
View(temp_pyr)
colnames(temp_pyr)
# Baseline of the HC_PYR thermogram
base_ch_pyr<-read.csv("base_ch_pyr.csv", header=T, row.names=1, sep=";", dec=",")
class(base_ch_pyr)
dim(base_ch_pyr)
View(base_ch_pyr)
base_ch_pyr
ch_pyr_corrected<-matrix(0,nrow(ch_pyr),ncol(ch_pyr))
for (i in 1:nrow(base_ch_pyr)) {
  ch_pyr_corrected[,i]<-ch_pyr[,i]-base_ch_pyr[i,1]
}
### Plot elution profiles
library(hyperSpec)
chpyr<-new("hyperSpec", spc = t(ch_pyr))
chpyr@wavelength
chpyr@label
chpyr@label$.wavelength<-"time"
chpyr@label$spc<-"I / a.u."
chpyr
chpyrcorrected<-new("hyperSpec", spc = t(ch_pyr_corrected))
chpyrcorrected@wavelength
chpyrcorrected@label
chpyrcorrected@label$.wavelength<-"time"
chpyrcorrected@label$spc<-"I / a.u."
chpyrcorrected

plot(chpyr, lines.args = list (pch = 20, lwd=2), wl.reverse= F, spc.nmax = ncol(ch_pyr), col="darkblue",  wl.range = c (min ~ max))
plot(chpyrcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, spc.nmax = ncol(ch_pyr), col="darkred",  wl.range = c (min ~ max), add=T)

library(pracma)
AreaCHpyr<-matrix(0,ncol(ch_pyr),90)
colnames(AreaCHpyr)<-c(1:90)
library(stringr)
colnames(AreaCHpyr)<-paste("t_ch_pyr_",colnames(AreaCHpyr),sep="")
rownames(AreaCHpyr)<-rownames(chpyr$spc)
View(AreaCHpyr)
for (i in 1:ncol(ch_pyr)) {
  for (j in 1:90) {
    AreaCHpyr[i,j]<-temp_pyr[200+which(cumtrapz(chpyrcorrected[,,200~nrow(ch_pyr)]@wavelength,matrix(chpyrcorrected[i,,200~nrow(ch_pyr)]$spc,nrow(ch_pyr)-199,1))>=trapz(chpyrcorrected[,,200~nrow(ch_pyr)]@wavelength,chpyrcorrected[i,,200~nrow(ch_pyr)]$spc)*j/100)[1],i]
  }}
View(AreaCHpyr)

###########################################################
### Calcul indices TLHC (Saenger et al., 2015) and I & R indices (Sebag et al. 2016) 
#   on the HC_PYR thermogram
library(pracma)
AreaCHpyr_IR<-matrix(0,ncol(ch_pyr),15)
rownames(AreaCHpyr_IR)<-rownames(chpyr$spc)
colnames(AreaCHpyr_IR)<-c("AREA_CH_PYR_total","Area_200_340_CH","Area_340_400_CH","Area_400_460_CH","Area_sup_460_CH","A1","A2","A3","A4","indice_i","indice_r","Area_inf_450_CH","Area_sup_450_CH","cl_ci","cp")
View(AreaCHpyr_IR)

for (i in 1:ncol(ch_pyr)) {
  AreaCHpyr_IR[i,1]<-trapz(chpyrcorrected[,,200~1260]@wavelength,chpyrcorrected[i,,200~1260]$spc)
  AreaCHpyr_IR[i,2]<-trapz(chpyrcorrected[,,200~which(temp_pyr[,i]>=340)[1]]@wavelength,chpyrcorrected[i,,200~which(temp_pyr[,i]>=340)[1]]$spc)
  AreaCHpyr_IR[i,3]<-trapz(chpyrcorrected[,,which(temp_pyr[,i]>=340)[1]~which(temp_pyr[,i]>=400)[1]]@wavelength,chpyrcorrected[i,,which(temp_pyr[,i]>=340)[1]~which(temp_pyr[,i]>=400)[1]]$spc)
  AreaCHpyr_IR[i,4]<-trapz(chpyrcorrected[,,which(temp_pyr[,i]>=400)[1]~which(temp_pyr[,i]>=460)[1]]@wavelength,chpyrcorrected[i,,which(temp_pyr[,i]>=400)[1]~which(temp_pyr[,i]>=460)[1]]$spc)
  AreaCHpyr_IR[i,5]<-trapz(chpyrcorrected[,,which(temp_pyr[,i]>=460)[1]~1260]@wavelength,chpyrcorrected[i,,which(temp_pyr[,i]>=460)[1]~1260]$spc)
  AreaCHpyr_IR[i,6]<-100*(AreaCHpyr_IR[i,2]/AreaCHpyr_IR[i,1])
  AreaCHpyr_IR[i,7]<-100*(AreaCHpyr_IR[i,3]/AreaCHpyr_IR[i,1])
  AreaCHpyr_IR[i,8]<-100*(AreaCHpyr_IR[i,4]/AreaCHpyr_IR[i,1])
  AreaCHpyr_IR[i,9]<-100*(AreaCHpyr_IR[i,5]/AreaCHpyr_IR[i,1])
  AreaCHpyr_IR[i,10]<-log10((AreaCHpyr_IR[i,6]+AreaCHpyr_IR[i,7])/AreaCHpyr_IR[i,8])
  AreaCHpyr_IR[i,11]<-(AreaCHpyr_IR[i,8]+AreaCHpyr_IR[i,9])/100
  AreaCHpyr_IR[i,12]<-trapz(chpyrcorrected[,,200~which(temp_pyr[,i]>=450)[1]]@wavelength,chpyrcorrected[i,,200~which(temp_pyr[,i]>=450)[1]]$spc)
  AreaCHpyr_IR[i,13]<-trapz(chpyrcorrected[,,which(temp_pyr[,i]>=450)[1]~1260]@wavelength,chpyrcorrected[i,,which(temp_pyr[,i]>=450)[1]~1260]$spc)
  AreaCHpyr_IR[i,14]<-AreaCHpyr_IR[i,12]/AreaCHpyr_IR[i,1]
  AreaCHpyr_IR[i,15]<-AreaCHpyr_IR[i,13]/AreaCHpyr_IR[i,1]
                  }
View(AreaCHpyr_IR)

###########################################################################################################
### Compute the temperature parameters for the CO2_PYR thermogram
############################################################################################################
rm(list=ls())
### Read files
# CO2_PYR thermogram
co2pyr<-read.csv("co2_pyr.csv", header=T, row.names=1, sep=";")
dim(co2pyr)
View(co2pyr)
#??? Oven temperature during pyrolysis
temp_pyr<-read.csv("temp_pyr.csv", header=T, row.names=1, sep=";")
dim(temp_pyr)
View(temp_pyr)
### Plot elution profiles
library(hyperSpec)
CO2pyr<-new("hyperSpec", spc = t(co2pyr))
CO2pyr@wavelength
CO2pyr@label
CO2pyr@label$.wavelength<-"time"
CO2pyr@label$spc<-"I / a.u."
CO2pyr
plot(CO2pyr, lines.args = list (pch = 20, lwd=2), wl.reverse= F, spc.nmax = ncol(co2pyr), col="darkblue", wl.range = c (min ~ max))

library(pracma)
AreaCO2pyr<-matrix(0,ncol(co2pyr),100)
colnames(AreaCO2pyr)<-c(1:100)
library(stringr)
colnames(AreaCO2pyr)<-paste("t_co2_pyr_",colnames(AreaCO2pyr),sep="")
rownames(AreaCO2pyr)<-rownames(CO2pyr$spc)
View(AreaCO2pyr)

for (i in 1:ncol(co2pyr)) {
  for (j in 1:100) {
    AreaCO2pyr[i,j]<-temp_pyr[200+which(cumtrapz(CO2pyr[,,200~which(temp_pyr[,i]>=560)[1]]@wavelength,matrix(CO2pyr[i,,200~which(temp_pyr[,i]>=560)[1]]$spc,which(temp_pyr[,i]>=560)[1]-199,1))>=trapz(CO2pyr[,,200~which(temp_pyr[,i]>=560)[1]]@wavelength,CO2pyr[i,,200~which(temp_pyr[,i]>=560)[1]]$spc)*j/100)[1],i]
  }}
View(AreaCO2pyr)

############################################################################################################
### Compute the temperature parameters for the CO_PYR thermogram
############################################################################################################
rm(list=ls())
### Read files
# CO_PYR thermogram
copyr<-read.csv("co_pyr.csv", header=T, row.names=1, sep=";")
dim(copyr)
View(copyr)
# Temperature of the oven during pyrolysis
temp_pyr<-read.csv("temp_pyr.csv", header=T, row.names=1, sep=";")
dim(temp_pyr)
View(temp_pyr)
### Plot elution profiles
library(hyperSpec)
COpyr<-new("hyperSpec", spc = t(copyr))
COpyr@wavelength
COpyr@label
COpyr@label$.wavelength<-"time"
COpyr@label$spc<-"I / a.u."
COpyr
plot(COpyr, lines.args = list (pch = 20, lwd=2), wl.reverse= F, spc.nmax = ncol(copyr), col="darkblue", wl.range = c (min ~ max))

library(pracma)
AreaCOpyr<-matrix(0,ncol(copyr),100)
colnames(AreaCOpyr)<-c(1:100)
library(stringr)
colnames(AreaCOpyr)<-paste("t_co_pyr_",colnames(AreaCOpyr),sep="")
rownames(AreaCOpyr)<-rownames(COpyr$spc)
View(AreaCOpyr)

for (i in 1:ncol(copyr)) {
  for (j in 1:100) {
    AreaCOpyr[i,j]<-temp_pyr[200+which(cumtrapz(COpyr[,,200~which(temp_pyr[,i]>=560)[1]]@wavelength,matrix(COpyr[i,,200~which(temp_pyr[,i]>=560)[1]]$spc,which(temp_pyr[,i]>=560)[1]-199,1))>=trapz(COpyr[,,200~which(temp_pyr[,i]>=560)[1]]@wavelength,COpyr[i,,200~which(temp_pyr[,i]>=560)[1]]$spc)*j/100)[1],i]
  }}
View(AreaCOpyr)

############################################################################################################
### Compute the temperature parameters for the CO2_OX thermogram
############################################################################################################
rm(list=ls())
### Read files
# CO2_OX thermogram
CO2_OX<-read.csv("co2_ox.csv", header=T, row.names=1, sep=";")
dim(CO2_OX)
View(CO2_OX)
colnames(CO2_OX)
# Temperature of the oven during oxidation
temp_ox<-read.csv("temp_ox.csv", header=T, row.names=1, sep=";")
dim(temp_ox)
View(temp_ox)
# Baseline of the CO2_OX thermogram
base_co2_ox<-read.csv("base_co2_ox.csv", header=T, row.names=1, sep=";", dec=",")
class(base_co2_ox)
dim(base_co2_ox)
View(base_co2_ox)
base_co2_ox
co2_ox_corrected<-matrix(0,nrow(CO2_OX),ncol(CO2_OX))
for (i in 1:nrow(base_co2_ox)) {
  co2_ox_corrected[,i]<-CO2_OX[,i]-base_co2_ox[i,1]
}
### Plot elution profiles
library(hyperSpec)
co2_OX<-new("hyperSpec", spc = t(CO2_OX))
co2_OX@wavelength
co2_OX@label
co2_OX@label$.wavelength<-"time"
co2_OX@label$spc<-"I / a.u."
co2_OX

co2_OXcorrected<-new("hyperSpec", spc = t(co2_ox_corrected))
co2_OXcorrected@wavelength
co2_OXcorrected@label
co2_OXcorrected@label$.wavelength<-"time"
co2_OXcorrected@label$spc<-"I / a.u."
co2_OXcorrected

plot(co2_OX, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(CO2_OX), wl.range = c (min ~ max))
plot(co2_OXcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(CO2_OX), wl.range = c (min ~ max), add=T)

library(pracma)
Areaco2OX<-matrix(0,ncol(CO2_OX),100)
colnames(Areaco2OX)<-c(1:100)
library(stringr)
colnames(Areaco2OX)<-paste("t_co2_ox_",colnames(Areaco2OX),sep="")
rownames(Areaco2OX)<-rownames(co2_OX$spc)
View(Areaco2OX)

for (i in 1:ncol(CO2_OX)) {
  for (j in 1:100) {
    Areaco2OX[i,j]<-temp_ox[which(cumtrapz(co2_OXcorrected[,,1~which(temp_ox[,i]>=611)[1]]@wavelength,matrix(co2_OXcorrected[i,,1~which(temp_ox[,i]>=611)[1]]$spc,which(temp_ox[,i]>=611)[1],1))>=trapz(co2_OXcorrected[,,1~which(temp_ox[,i]>=611)[1]]@wavelength,co2_OXcorrected[i,,1~which(temp_ox[,i]>=611)[1]]$spc)*j/100)[1],i]
  }}
View(Areaco2OX)

############################################################################################################
### Compute the temperature parameters for the CO_OX thermogram
############################################################################################################
rm(list=ls())
### Read files
# CO_OX thermogram
coOX<-read.csv("co_ox.csv", header=T, row.names=1, sep=";")
dim(coOX)
View(coOX)
# Temperature of the oven during oxidation
temp_ox<-read.csv("temp_ox.csv", header=T, row.names=1, sep=";")
dim(temp_ox)
View(temp_ox)
# Baseline of the CO_OX thermogram
base_co_ox<-read.csv("base_co_ox.csv", header=T, row.names=1, sep=";", dec=",")
class(base_co_ox)
dim(base_co_ox)
View(base_co_ox)
base_co_ox
co_ox_corrected<-matrix(0,nrow(coOX),ncol(coOX))

for (i in 1:nrow(base_co_ox)){
  co_ox_corrected[,i]<-coOX[,i]-base_co_ox[i,1]
}

### Plot elution profiles
library(hyperSpec)
COOX<-new("hyperSpec", spc = t(coOX))
COOX@wavelength
COOX@label
COOX@label$.wavelength<-"time"
COOX@label$spc<-"I / a.u."
COOX
COOXcorrected<-new("hyperSpec", spc = t(co_ox_corrected))
COOXcorrected@wavelength
COOXcorrected@label
COOXcorrected@label$.wavelength<-"time"
COOXcorrected@label$spc<-"I / a.u."
COOXcorrected
plot(COOX, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(coOX), wl.range = c (min ~ max))
plot(COOXcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(coOX), wl.range = c (min ~ max), add=T)

library(pracma)
AreacoOX<-matrix(0,ncol(coOX),90)
colnames(AreacoOX)<-c(1:90)
library(stringr)
colnames(AreacoOX)<-paste("t_co_ox_",colnames(AreacoOX),sep="")
rownames(AreacoOX)<-rownames(COOX$spc)
View(AreacoOX)

for (i in 1:ncol(coOX)) {
  for (j in 1:90) {
    AreacoOX[i,j]<-temp_ox[which(cumtrapz(COOXcorrected[,,0~nrow(coOX)]@wavelength,matrix(COOXcorrected[i,,0~nrow(coOX)]$spc,nrow(coOX),1))>=trapz(COOXcorrected[,,0~nrow(coOX)]@wavelength,COOXcorrected[i,,0~nrow(coOX)]$spc)*j/100)[1],i]
  }}
View(AreacoOX)

##############################################################################
### Computing PseudoS1 and S2 parameters 
##############################################################################
### Pre-processing the signal before that...
### Correction "carbon-equivalent" of thermograms
# read files
# kFID parameter (used to adapt CH_pyr thermogram)
kfid<-read.csv("kFID.csv", header=T, row.names=1, sep=";", dec=",")
class(kfid)
dim(kfid)
View(kfid)
kfid
# Sample weight used to normalize all thermograms
sample_weight<-read.csv("sample_weight.csv", header=T, row.names=1, sep=";", dec=",")
class(sample_weight)
dim(sample_weight)
View(sample_weight)
sample_weight

### Normalizing HC_PYR thermogram (baseline corrected)
# Raw thermogram is 
# multiplied by kfid, by 100, by the relative molar mass (0.83) 
# divided by sample weight
ch_pyr_corrected_carbon<-matrix(0,nrow(ch_pyr_corrected),ncol(ch_pyr_corrected))
for (i in 1:nrow(base_ch_pyr)){
  ch_pyr_corrected_carbon[,i]<-ch_pyr_corrected[,i]*kfid[i,1]*100*0.83/sample_weight[i,1]
}

chpyrcorrectedcarbon<-new("hyperSpec", spc = t(ch_pyr_corrected_carbon))
chpyrcorrectedcarbon@wavelength
chpyrcorrectedcarbon@label
chpyrcorrectedcarbon@label$.wavelength<-"time"
chpyrcorrectedcarbon@label$spc<-"I / a.u."
chpyrcorrectedcarbon

plot(chpyrcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))
plot(chpyrcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))

write.csv2(chpyrcorrectedcarbon,"CHPYR_eqC.csv")

### AREA S1 & S2
library(pracma)
AreaS1pyr<-matrix(0,ncol(temp_pyr),1)
colnames(AreaS1pyr)<-c("S1area")
rownames(AreaS1pyr)<-colnames(temp_pyr)
View(AreaS1pyr)

for (i in 1:ncol(temp_pyr)) {
  AreaS1pyr[i,1]<-trapz(chpyrcorrectedcarbon[,,0~200]@wavelength,chpyrcorrectedcarbon[i,,0~200]$spc)
}

library(pracma)
AreaS2pyr<-matrix(0,ncol(temp_pyr),1)
colnames(AreaS2pyr)<-c("S2area")
rownames(AreaS2pyr)<-colnames(temp_pyr)
View(AreaS2pyr)

for (i in 1:ncol(temp_pyr)) {
  AreaS2pyr[i,1]<-trapz(chpyrcorrectedcarbon[,,200~nrow(temp_pyr)]@wavelength,chpyrcorrectedcarbon[i,,200~nrow(temp_pyr)]$spc)
}
View(AreaS2pyr)
length(AreaS2pyr)

### Normalizing CO_PYR thermogram
# Thermogram is multiplied the relative molar mass of C 12/28,
# divided by 1000,
# divided by sample weight
co_pyr_carbon<-matrix(0,nrow(copyr),ncol(copyr))
for (i in 1:ncol(temp_pyr)){
  co_pyr_carbon[,i]<-copyr[,i]*(12/28)/1000/sample_weight[i,1]
}

copyrcorrectedcarbon<-new("hyperSpec", spc = t(co_pyr_carbon))
copyrcorrectedcarbon@wavelength
copyrcorrectedcarbon@label
copyrcorrectedcarbon@label$.wavelength<-"time"
copyrcorrectedcarbon@label$spc<-"I / a.u."
copyrcorrectedcarbon

plot(COpyr, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))
plot(copyrcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))

### S3CO
library(pracma)
AreaS3copyr<-matrix(0,ncol(temp_pyr),1)
colnames(AreaS3copyr)<-c("S3coarea")
rownames(AreaS3copyr)<-colnames(temp_pyr)
View(AreaS3copyr)

for (i in 1:ncol(temp_pyr)) {
  AreaS3copyr[i,1]<-trapz(copyrcorrectedcarbon[,,0~nrow(temp_pyr)]@wavelength,copyrcorrectedcarbon[i,,0~nrow(temp_pyr)]$spc)
}
View(AreaS3copyr)
length(AreaS3copyr)

### Normalizing CO2_PYR thermogram
# Thermogram is multiplied by the relative molar mass C 12/44,
# divided by 1000,
# divided by sample weight

co2_pyr_carbon<-matrix(0,nrow(co2pyr),ncol(co2pyr))
for (i in 1:ncol(temp_pyr)) {
  co2_pyr_carbon[,i]<-co2pyr[,i]*(12/44)/1000/sample_weight[i,1]
}

co2pyrcorrectedcarbon<-new("hyperSpec", spc = t(co2_pyr_carbon))
co2pyrcorrectedcarbon@wavelength
co2pyrcorrectedcarbon@label
co2pyrcorrectedcarbon@label$.wavelength<-"time"
co2pyrcorrectedcarbon@label$spc<-"I / a.u."
co2pyrcorrectedcarbon

plot(CO2pyr, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))
plot(co2pyrcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(ch_pyr_corrected_carbon), wl.range = c (min ~ max))

### S3CO2
library(pracma)
AreaS3co2pyr<-matrix(0,ncol(temp_pyr),1)
colnames(AreaS3co2pyr)<-c("S3co2area")
rownames(AreaS3co2pyr)<-colnames(temp_pyr)
View(AreaS3co2pyr)

for (i in 1:ncol(temp_pyr)) {
  AreaS3co2pyr[i,1]<-trapz(co2pyrcorrectedcarbon[,,0~nrow(temp_pyr)]@wavelength,co2pyrcorrectedcarbon[i,,0~nrow(temp_pyr)]$spc)
}
View(AreaS3co2pyr)
length(AreaS3co2pyr)

### Normalizing CO_OX thermogram (baseline corrected)
# Thermogram is multiplied by relative molar mass C 12/28,
# divided by 1000,
# divided by sample weight

co_ox_corrected_carbon<-matrix(0,nrow(co_ox_corrected),ncol(co_ox_corrected))
for (i in 1:nrow(base_co_ox))  {
  co_ox_corrected_carbon[,i]<-co_ox_corrected[,i]*(12/28)/1000/sample_weight[i,1]
}

cooxcorrectedcarbon<-new("hyperSpec", spc = t(co_ox_corrected_carbon))
cooxcorrectedcarbon@wavelength
cooxcorrectedcarbon@label
cooxcorrectedcarbon@label$.wavelength<-"time"
cooxcorrectedcarbon@label$spc<-"I / a.u."
cooxcorrectedcarbon

plot(COOXcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(co_ox_corrected_carbon), wl.range = c (min ~ max))
plot(cooxcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(co_ox_corrected_carbon), wl.range = c (min ~ max))

### S4CO
library(pracma)
AreaS4coox<-matrix(0,ncol(temp_ox),1)
colnames(AreaS4coox)<-c("S4coarea")
rownames(AreaS4coox)<-colnames(temp_ox)
View(AreaS4coox)

for (i in 1:ncol(temp_ox)) {
  AreaS4coox[i,1]<-trapz(cooxcorrectedcarbon[,,0~nrow(temp_ox)]@wavelength,cooxcorrectedcarbon[i,,0~nrow(temp_ox)]$spc)
}
View(AreaS4coox)
length(AreaS4coox)

### Normalizing CO2_OX thermogram (baseline corrected)
# Thermogram is multiplied by relative molar mass C 12/44,
# divided by 1000,
# divided by sample weight

co2_ox_corrected_carbon<-matrix(0,nrow(co2_ox_corrected),ncol(co2_ox_corrected))
for (i in 1:nrow(base_co2_ox)) {
  co2_ox_corrected_carbon[,i]<-co2_ox_corrected[,i]*(12/44)/1000/sample_weight[i,1]
}

co2oxcorrectedcarbon<-new("hyperSpec", spc = t(co2_ox_corrected_carbon))
co2oxcorrectedcarbon@wavelength
co2oxcorrectedcarbon@label
co2oxcorrectedcarbon@label$.wavelength<-"time"
co2oxcorrectedcarbon@label$spc<-"I / a.u."
co2oxcorrectedcarbon

plot(co2_OXcorrected, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = ncol(co2_ox_corrected_carbon), wl.range = c (min ~ max))
plot(co2oxcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(co2_ox_corrected_carbon), wl.range = c (min ~ max))
spc.identify(plot(co2oxcorrectedcarbon, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkred", spc.nmax = ncol(co2_ox_corrected_carbon), wl.range = c (min ~ max)))

### S5CO2
library(pracma)
AreaS5co2ox<-matrix(0,ncol(temp_ox),1)
colnames(AreaS5co2ox)<-c("S5co2area")
rownames(AreaS5co2ox)<-colnames(temp_ox)
View(AreaS5co2ox)

for (i in 1:ncol(temp_ox)) {
  AreaS5co2ox[i,1]<-trapz(co2oxcorrectedcarbon[,,0~nrow(temp_ox)]@wavelength,co2oxcorrectedcarbon[i,,0~nrow(temp_ox)]$spc)
}
View(AreaS5co2ox)
length(AreaS5co2ox)

###################################################
### Summing the thermograms of the pyrolysis phase (SUMPYR)
### (=HC_PYR+CO_PYR+CO2_PYR in "carbon-equivalent" signals; mgC/g/s)
### SUMPYR
SUMPYR<-chpyrcorrectedcarbon+copyrcorrectedcarbon+co2pyrcorrectedcarbon
SUMPYR
SUMPYR@label$.wavelength<-""
SUMPYR@label$spc<-""

rownames(SUMPYR)<-colnames(temp_pyr)
write.csv2(SUMPYR,"SUMPYR_eqC.csv")

plot(SUMPYR[,,0~200], lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = 850, wl.range = c (min ~ max))
plot(SUMPYR[,,200~1300], lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = 850, wl.range = c (min ~ max))

### Calculate instantaneous proportion of C evoluted as HC during pyrolysis
### Dividing "CHpyrcorrectedcarbon" by "SUMPYR"
prop_hc<-matrix(0,nrow(chpyrcorrectedcarbon$spc),ncol(chpyrcorrectedcarbon$spc))
dim(prop_hc)
colnames(prop_hc)<-colnames(chpyr$spc)
rownames(prop_hc)<-rownames(chpyr$spc)
View(prop_hc)

for(i in 1:nrow(chpyrcorrectedcarbon$spc)){
prop_hc[i,]<-sweep(chpyrcorrectedcarbon[i,,], 2, SUMPYR[i,,], "/")$spc
}

View(prop_hc)
write.csv2(prop_hc,"proportion_HC_pyr_instant.csv")

#######################
### Summing the thermograms of the oxidation phase (SUMOX)
### (=CO_OX+CO2OX in "carbon-equivalent" signals; mgC/g/s
SUMOX<-cooxcorrectedcarbon+co2oxcorrectedcarbon
SUMOX
SUMOX@label$.wavelength<-"Time"
SUMOX@label$spc<-"mgC/gsoil"
rownames(SUMOX)<-colnames(temp_ox)
write.csv2(SUMOX,"SUMOX_eqC.csv")

spc.identify(plot(SUMOX, lines.args = list (pch = 20, lwd=2), wl.reverse= F, col="darkblue", spc.nmax = 850, wl.range = c (min ~ max)))
rownames(SUMOX)

### File with areas and ratios
# To compute CTOT en mgC/g, two methods are compared: 
# all signal PYR+OX or PYR+OX cut at 611°C (to reduce noise for samples low in C)
library(pracma)
AreaSUM<-matrix(0,ncol(temp_pyr),15)
colnames(AreaSUM)<-c("Area_PYR_actual_S1_mgC_g","Area_PYR_actual_S2_mgC_g",
                     "AreaSUMPYR_pseudoS1_mgC_g","AreaSUMPYR_pseudoS2_mgC_g",
                     "AreaSUMPYR_total_mgC_g","AreaSUMOX_mgC_g",
                     "AreaSUMOX_mgC_g_cutox611","AreaTOTAL_mgC_g",
                     "AreaTOTAL_mgC_g_cutox611","ratio_pseudoS1_Ctot",
                     "ratio_Cpyr_Ctot","ratio_pseudoS1_Ctot_cutox611",
                     "ratio_Cpyr_Ctot_cutox611","ratio_S2_Cpyr",
                     "ratio_S1S2_Cpyr")
rownames(AreaSUM)<-colnames(temp_pyr)
View(AreaSUM)

for (i in 1:ncol(temp_pyr)) {
  AreaSUM[i,1]<-AreaS1pyr[i]
  AreaSUM[i,2]<-AreaS2pyr[i]
  AreaSUM[i,3]<-trapz(SUMPYR[,,0~200]@wavelength,SUMPYR[i,,0~200]$spc)
  AreaSUM[i,4]<-trapz(SUMPYR[,,200~nrow(temp_pyr)]@wavelength,SUMPYR[i,,200~nrow(temp_pyr)]$spc)
  AreaSUM[i,5]<-AreaSUM[i,3]+AreaSUM[i,4]
  AreaSUM[i,6]<-trapz(SUMOX[,,0~nrow(temp_ox)]@wavelength,SUMOX[i,,0~nrow(temp_ox)]$spc)
  AreaSUM[i,7]<-trapz(SUMOX[,,0~which(temp_ox[,i]>=611)[1]]@wavelength,SUMOX[i,,0~which(temp_ox[,i]>=611)[1]]$spc)
  AreaSUM[i,8]<-AreaSUM[i,5]+AreaSUM[i,6]
  AreaSUM[i,9]<-AreaSUM[i,5]+AreaSUM[i,7]
  AreaSUM[i,10]<-AreaSUM[i,3]/AreaSUM[i,8]
  AreaSUM[i,11]<-AreaSUM[i,5]/AreaSUM[i,8]
  AreaSUM[i,12]<-AreaSUM[i,3]/AreaSUM[i,9]
  AreaSUM[i,13]<-AreaSUM[i,5]/AreaSUM[i,9]
  AreaSUM[i,14]<-AreaSUM[i,2]/AreaSUM[i,5]
  AreaSUM[i,15]<-(AreaSUM[i,1]+AreaSUM[i,2])/AreaSUM[i,5]
}

View(AreaSUM)
write.csv2(AreaSUM,"AREAS_eqCARBON.csv")

###################################################################
### Calculating temperature parameters 
### SUMPYR : version starting after the isotherm (200s) at 200°C

library(pracma)
AreaPYR<-matrix(0,ncol(temp_pyr),90)
colnames(AreaPYR)<-c(1:90)
library(stringr)
colnames(AreaPYR)<-paste("t_PYR_",colnames(AreaPYR),sep="")
rownames(AreaPYR)<-colnames(temp_pyr)
View(AreaPYR)

for (i in 1:ncol(temp_pyr)) {
  for (j in 1:90) {
    AreaPYR[i,j]<-temp_pyr[200+which(cumtrapz(SUMPYR[,,200~nrow(temp_pyr)]@wavelength,matrix(SUMPYR[i,,200~nrow(temp_pyr)]$spc,nrow(temp_pyr)-199,1))>=trapz(SUMPYR[,,200~nrow(temp_pyr)]@wavelength,SUMPYR[i,,200~nrow(temp_pyr)]$spc)*j/100)[1],i]
  }}

View(AreaPYR)
write.csv2(AreaPYR,"t_parameters_SUMPYR_200s_to_full_signal.csv")

####################################################################
### Calculating temperature parameters 
### SUMPYR : version starting at the onset of pyrolysis

library(pracma)
AreaPYR2<-matrix(0,ncol(temp_pyr),100)
colnames(AreaPYR2)<-c(1:100)
library(stringr)
colnames(AreaPYR2)<-paste("t_pyr_full_signal_",colnames(AreaPYR2),sep="")
rownames(AreaPYR2)<-colnames(temp_pyr)
View(AreaPYR2)
for (i in 1:ncol(temp_pyr)) {
  for (j in 1:100) {
    AreaPYR2[i,j]<-temp_pyr[which(cumtrapz(SUMPYR[,,0~nrow(temp_pyr)]@wavelength,matrix(SUMPYR[i,,0~nrow(temp_pyr)]$spc,nrow(temp_pyr),1))>=trapz(SUMPYR[,,0~nrow(temp_pyr)]@wavelength,SUMPYR[i,,0~nrow(temp_pyr)]$spc)*j/100)[1],i]
  }}
View(AreaPYR2)

### ATTENTION : 
# Replace values above 645 degrees by 650 
for (i in 1:ncol(temp_pyr)) {
  AreaPYR2[i,c(which(AreaPYR2[i,]>=645)[1]:100)]<-650
}
View(AreaPYR2)
write.csv2(AreaPYR2,"t_parameters_SUMPYR_full_signal.csv")

par(mar=c(5,5,4,5)+.1)
plot(seq(200,650,length.out=100),c(1:100), type="n", las=1, ylab="Pyrolyzed C (%)", xlab="T (°C)", cex.lab=1.5, cex.axis=1.25, xlim=c(200,650), ylim=c(0,100))
for (i in 1:nrow(AreaPYR2)) {
  lines(AreaPYR2[i,c(1:100)],c(1:100), lwd=4, type="l")
}

####################################################################
### PLOTTING proportion de C évolué sous forme hydrocarbone durant la phase de pyrolyse 
### en proportion of total pyrolyzed C
library(pracma)
samp<-matrix(0,ncol(temp_pyr),100)
temp<-matrix(0,ncol(temp_pyr),100)
AreaPYR3<-matrix(0,ncol(temp_pyr),100)
colnames(samp)<-c(1:100)
colnames(temp)<-c(1:100)
colnames(AreaPYR3)<-c(1:100)
library(stringr)
colnames(AreaPYR3)<-paste("prop_HC_cumul_",colnames(AreaPYR3),sep="")
rownames(AreaPYR3)<-colnames(temp_pyr)
View(AreaPYR3)

for (i in 1:ncol(temp_pyr)) {
samp[i,1]<-2
temp[i,1]<-200
for (j in 1:99)
{
  samp[i,j+1]<-which(temp_pyr[,i]>=seq(temp_pyr[1,i], temp_pyr[which(temp_pyr[,i]>=649)[1],i], length.out = 100)[j+1])[1]
  temp[i,j+1]<-temp_pyr[which(temp_pyr[,i]>=seq(temp_pyr[1,i], temp_pyr[which(temp_pyr[,i]>=649)[1],i], length.out = 100)[j+1])[1],i]
}}

samp
temp

for (i in 1:ncol(temp_pyr))  {
    test<-cumtrapz(chpyrcorrectedcarbon[,,0~nrow(temp_pyr)]@wavelength,matrix(chpyrcorrectedcarbon[i,,0~nrow(temp_pyr)]$spc,nrow(temp_pyr),1))/
    cumtrapz(SUMPYR[,,0~nrow(temp_pyr)]@wavelength,matrix(SUMPYR[i,,0~nrow(temp_pyr)]$spc,nrow(temp_pyr),1))
    AreaPYR3[i,]<-test[c(samp[i,]),]
  }

View(AreaPYR3)
colnames(AreaPYR3)<-temp[1,]
View(AreaPYR3)
write.csv2(AreaPYR3,"proportion_HC_pyr_cumul.csv")

par(new=TRUE)
for (i in 1:nrow(AreaPYR3)) {
  lines(temp[i,], AreaPYR3[i,]*100, lwd=3, lty=3, type="l", xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0,100))
}
axis(4, las=1, cex.lab=1.5, cex.axis=1.25)
mtext("Hydrocarbons (%)", side=4, line=3, cex = 1.5)
legend("bottomright", col=c("black"), lwd=c(4,3), lty=c(1,3), cex=1.5, legend=c("Pyrolyzed C","Hydrocarbons"))

####################################################################
### PLOTTING evolved C during pyrolysis (% of total C) = f (T)
library(pracma)
AreaPYR4<-matrix(0,ncol(temp_pyr),100)
colnames(AreaPYR4)<-c(1:100)
library(stringr)
colnames(AreaPYR4)<-paste("t_pyr_fulsign_totC_",colnames(AreaPYR4),sep="")
rownames(AreaPYR4)<-colnames(temp_pyr)
View(AreaPYR4)
for (i in 1:ncol(temp_pyr)) {
  for (j in 1:100) {
    AreaPYR4[i,j]<-temp_pyr[which(cumtrapz(SUMPYR[,,0~nrow(temp_pyr)]@wavelength,matrix(SUMPYR[i,,0~nrow(temp_pyr)]$spc,nrow(temp_pyr),1))>=AreaSUM[i,9]*j/100)[1],i]
  }}
View(AreaPYR4)
write.csv2(AreaPYR4,"t_parameters_SUMPYR_fulsig_totC.csv")

# Tracing plot till 629°C 
par(mar=c(5,5,4,5)+.1)
plot(seq(200,650,length.out=100),c(1:100), type="n", las=1, ylab="Pyrolyzed C (% of total C)", xlab="T (°C)", cex.lab=1.5, cex.axis=1.25, xlim=c(200,650), ylim=c(0,100))
for (i in 1:nrow(AreaPYR4)) {
  lines(AreaPYR4[i,c(1:which(AreaPYR4[i,]>=629)[1])],c(1:which(AreaPYR4[i,]>=629)[1]), lwd=4, type="l")
}

par(new=TRUE)
for (i in 1:nrow(AreaPYR3)) {
  lines(temp[i,], AreaPYR3[i,]*100, lwd=3, lty=3, type="l", xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0,100))
}
axis(4, las=1, cex.lab=1.5, cex.axis=1.25)
mtext("Hydrocarbons (%)", side=4, line=3, cex = 1.5)
legend("topright", col=c("black"), lwd=c(4,3), lty=c(1,3), cex=1.5, legend=c("Pyrolyzed C","Hydrocarbons"))

####################################################################
### SUMOX
library(pracma)
AreaOX<-matrix(0,ncol(temp_ox),90)
colnames(AreaOX)<-c(1:90)
library(stringr)
colnames(AreaOX)<-paste("t_OX_",colnames(AreaOX),sep="")
rownames(AreaOX)<-colnames(temp_ox)
View(AreaOX)

for (i in 1:ncol(temp_ox)) {
  for (j in 1:90) {
    AreaOX[i,j]<-temp_ox[which(cumtrapz(SUMOX[,,0~which(temp_ox[,i]>=611)[1]]@wavelength,matrix(SUMOX[i,,0~which(temp_ox[,i]>=611)[1]]$spc,which(temp_ox[,i]>=611)[1],1))>=trapz(SUMOX[,,0~which(temp_ox[,i]>=611)[1]]@wavelength,SUMOX[i,,0~which(temp_ox[,i]>=611)[1]]$spc)*j/100)[1],i]
  }}

View(AreaOX)
write.csv2(AreaOX,"t_parameters_SUMOX_cutox611.csv")
