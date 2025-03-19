# Script for evaluating SDMs selecting pseudo-absences within a 20-kilometer disk arround the species observations

# Clear the environment
rm(list=ls())
# Load necessary libraries
library(raster)
library(tidyverse)
library(sp)
library(ecospat)
library(Metrics)

setwd("C:/Users/mofal/Desktop/R-Publication/Data")

# Load species names
species <- read.csv("Speciesnames.csv")
species <- species$species
spec1 <- species[1]

# Create data frame to store evaluation results
metric_final <- data.frame(value = as.numeric(), metric = as.character(), scale = as.character(), species = as.character(), model = as.character())

# Loop through species
for (spec1 in species){

# Change species names to fit data files
spec <- gsub("-",".",spec1, fixed=TRUE)

# Prepare file to load species observations
input_file_alps <- paste("moritz_johannes_20240416/disk/data_model/",spec,".csv",sep="")
input_file_mountains <- paste("Input_data_all_alps/",spec1,"_all",".csv",sep="")

# Load species observations at extent of the Alps
input_data_alps <- read.csv(input_file_alps, sep=",")
# Select species observations which were not used during model calibration
input_data_alps <- filter(input_data_alps, param == 0)

# Load species observations at extent of the Mountains
input_data_mountains <- read.csv(input_file_mountains, sep=";")
# Transform the coordinates to numeric
input_data_mountains$X_round <- as.numeric(gsub(",",".",input_data_mountains$X_round, fixed=TRUE))
input_data_mountains$Y_round <- as.numeric(gsub(",",".",input_data_mountains$Y_round, fixed=TRUE))

# Create paths to load binary and continuous model predictions
disk20_binary <- paste("C:/Users/mofal/Desktop/R-Publication/Data/moritz_johannes_20240416/disk/",spec,"/proj_current_","/proj_current__", spec,"_ensemble_TSSbin.tif",sep="")
disk20_continuous <- paste("C:/Users/mofal/Desktop/R-Publication/Data/moritz_johannes_20240416/disk/",spec,"/proj_current_","/proj_current__", spec,"_ensemble.tif",sep="")

# Load binary and continuous model predictions
prediction_binary <- raster::raster(disk20_binary) 
prediction_continuous <- raster(disk20_continuous)/1000 

# Select coordinates of presence observations at spatial extent of Alps
alps <- input_data_alps %>%dplyr::select(X_round, Y_round)

# Select presence observations at spatial extent of Schrankogel
S <- input_data_mountains %>% filter( mountains == "S") %>%
  dplyr::select(X_round, Y_round)

# Select presence observations at spatial extent of Racherin
G <- input_data_mountains %>% filter( mountains == "G") %>%
  dplyr::select(X_round, Y_round)

# Select presence observations at spatial extent of Hochschwab
H <- input_data_mountains %>% filter( mountains == "H") %>%
  dplyr::select(X_round, Y_round)

### Sample new pseudo-absences for model evaluation
# Extract cells for which we have information about the species presence
cells_1 <- data.frame(raster::extract(prediction_binary, input_data_mountains[1:2], cellnumber=T))
# Extract cell numbers
cells_1 <- cells_1[1]

# Create new raster layer which fits the binary predictions
pot_absences <- prediction_binary
for(i in 1:nrow(cells_1)){
  # set all cells where we have information about the species presence to NA
  pot_absences[cells_1[i,]] <- NA 
  # = Raster, where we have NA if species is present
  # later we will sample pseudo-absences in all cells !NA
}

# Lets evaluate!  
# True positive (TP), False positive (FP), False negatives (FN)
# True positive (TP): Data has presences and model predicts presence
# False negative (FN): Data has presence and model predicts absence
# False positive (FP): Data has absence and model predicts presence
# True negative (TN): Data has absence and model predicts absence
#
# Sensitivity (SN): TP/(TP+FN)
# Specificity (SP): TN/(TN+FP)
# True skill statistic (TSS): Sn + Sp - 1

### 1) Evaluation at spatial extent of the Alps
# Extract binary predictions at extent of the Alps
pres_alps <- data.frame(predicted= raster::extract(prediction_binary, alps))
pres_alps <- data.frame(observed = 1, predicted =pres_alps$predicted)

# Extract continuous predictions at extent of the Alps
pres_alps_cont <- data.frame(predicted= raster::extract(prediction_continuous, alps))
pres_alps_cont <- data.frame(observed = 1, predicted =pres_alps_cont$predicted)

# Extract True positive and False negative binary predictions
TP_alps <- sum(pres_alps$predicted == 1)/nrow(pres_alps)
FN_alps <- sum(pres_alps$predicted == 0)/nrow(pres_alps)

# Prepare dataframes to store evaluation according to TSS and AUC
TSS_alps_5 <- data.frame(value = as.numeric())
AUC_alps_5 <- data.frame(value = as.numeric())
# Evaluation of TSS and AUC with five iterations
for (t in 1:5){
  # Sample random pseudo-absences within the layer created before
  absent <- data.frame(x=numeric(), y=numeric(), layer=numeric())
  # Sample as many pseudo-absences as species obserations used for model calibration
  for (i in 1:nrow(alps)){
    # Create a 20 km buffer around each presence observation
    point <- alps %>% slice(i) 
    coordinates(point) <- ~X_round + Y_round
    crs(point) <- crs(prediction_binary)
    point <- spTransform(point, crs(prediction_binary))
    point <- raster::buffer(point, width=20000)
    
    # Crop buffer to raster which excludes presences
    pot_ab <- crop(pot_absences, point)
    pot_ab <- mask(pot_ab, point)
    
    # Sample pseudo-absences where there are no presence observation (= !NA)
    absence <- data.frame(sampleRandom(pot_ab, 1, na.rm=T, cells=F, xy=T))
    absent <- rbind(absent,absence)
  }
  # Extract coordinates of sampled pseudo-absences
  absences<- absent[1:2]
  
  # Extract binary predictions of the pseudo-absences
  abs_alps <- data.frame(predicted = raster::extract(prediction_binary, absences))
  abs_alps <- data.frame(observed = 0, predicted =abs_alps$predicted)
  
  # Extract binary predictions of the pseudo-absences
  abs_alps_cont <- data.frame(predicted = raster::extract(prediction_continuous, absences))
  abs_alps_cont <- data.frame(observed = 0, predicted =abs_alps_cont$predicted)
  
  # Extract False positive and True negative binary predictions
  FP_alps <- sum(abs_alps$predicted == 1)/nrow(abs_alps)
  TN_alps <- sum(abs_alps$predicted == 0)/nrow(abs_alps) 
  
  # Create data frame with observed vs predicted values (Binary)
  obs_pred_alps <- rbind(pres_alps,abs_alps)
  
  # Create data frame with observed vs predicted values (Continuous)
  obs_pred_alps_cont <- rbind(pres_alps_cont,abs_alps_cont)
  
  
  # Calculate Sensitivity, Specificity, True skill statistic and Area under the curve at extent of the Alps
  SENSI_alps <- data.frame(value = TP_alps/(TP_alps+FN_alps))
  SPECI_alps <- data.frame(value = TN_alps/(TN_alps+FP_alps))
  TSS_alps <- data.frame(value = SENSI_alps[1] + SPECI_alps[1] - 1)
  AUC_alps <- data.frame(value = auc(obs_pred_alps_cont$observed, obs_pred_alps_cont$predicted))
  
  # stack the results of TSS and AUC for the five evaluation runs!
  TSS_alps_5 <- rbind(TSS_alps_5,TSS_alps)
  AUC_alps_5 <- rbind(AUC_alps_5,AUC_alps)
}

# Create data frame with the results of TSS evaluation (=Mean of 5 runs)
TSS_alps <-  data.frame(value = mean(TSS_alps_5$value), metric= "TSS", scale = "ALPS", species = spec1)

# Create data frame with the results of AUC evaluation (=Mean of 5 runs)
AUC_alps <- data.frame(value = mean(AUC_alps_5$value), metric= "AUC", scale = "ALPS", species = spec1)

# Evaluation with the Boyce-index an presence only evaluation metric (1 run)
boyce_alps <- ecospat.boyce(prediction_continuous, alps, nclass = 0, window.w="default", res=1000, PEplot = FALSE,method = 'spearman')
boyce_alps <- boyce_alps$cor
# Create data frame with the results of Boyce-index
boyce_alps <- data.frame(value = boyce_alps, metric= "BOYCE", scale = "ALPS", species = spec1)


# Create data frame with the the results of TSS, AUC and Boyce-index
metric_alps <- rbind(TSS_alps,AUC_alps,boyce_alps)

### 2) Evaluation at spatial extent of the Schrankogel
# Extract binary predictions at extent of the S
pres_S <- data.frame(predicted= raster::extract(prediction_binary, S))
pres_S <- data.frame(observed = 1, predicted =pres_S$predicted)

# Extract continuous predictions at extent of the S
pres_S_cont <- data.frame(predicted= raster::extract(prediction_continuous, S))
pres_S_cont <- data.frame(observed = 1, predicted =pres_S_cont$predicted)

# Extract True positive and False negative binary predictions
TP_S <- sum(pres_S$predicted == 1)/nrow(pres_S)
FN_S <- sum(pres_S$predicted == 0)/nrow(pres_S)

# Prepare dataframes to store evaluation according to TSS and AUC
TSS_S_5 <- data.frame(value = as.numeric())
AUC_S_5 <- data.frame(value = as.numeric())
# Evaluation of TSS and AUC with five iterations
for (t in 1:5){
  # Sample random pseudo-absences within the layer created before
  absent <- data.frame(x=numeric(), y=numeric(), layer=numeric())
  # Sample as many pseudo-absences as species obserations used for model calibration
  for (i in 1:nrow(S)){
    # Create a 20 km buffer around each presence observation
    point <- S %>% slice(i) 
    coordinates(point) <- ~X_round + Y_round
    crs(point) <- crs(prediction_binary)
    point <- spTransform(point, crs(prediction_binary))
    point <- raster::buffer(point, width=20000)
    
    # Crop buffer to raster which excludes presences
    pot_ab <- crop(pot_absences, point)
    pot_ab <- mask(pot_ab, point)
    
    # Sample pseudo-absences where there are no presence observation (= !NA)
    absence <- data.frame(sampleRandom(pot_ab, 1, na.rm=T, cells=F, xy=T))
    absent <- rbind(absent,absence)
  }
  # Extract coordinates of sampled pseudo-absences
  absences<- absent[1:2]
  
  # Extract binary predictions of the pseudo-absences
  abs_S <- data.frame(predicted = raster::extract(prediction_binary, absences))
  abs_S <- data.frame(observed = 0, predicted =abs_S$predicted)
  
  # Extract binary predictions of the pseudo-absences
  abs_S_cont <- data.frame(predicted = raster::extract(prediction_continuous, absences))
  abs_S_cont <- data.frame(observed = 0, predicted =abs_S_cont$predicted)
  
  # Extract False positive and True negative binary predictions
  FP_S <- sum(abs_S$predicted == 1)/nrow(abs_S)
  TN_S <- sum(abs_S$predicted == 0)/nrow(abs_S) 
  
  # Create data frame with observed vs predicted values (Binary)
  obs_pred_S <- rbind(pres_S,abs_S)
  
  # Create data frame with observed vs predicted values (Continuous)
  obs_pred_S_cont <- rbind(pres_S_cont,abs_S_cont)
  
  
  # Calculate Sensitivity, Specificity, True skill statistic and Area under the curve at extent of the S
  SENSI_S <- data.frame(value = TP_S/(TP_S+FN_S))
  SPECI_S <- data.frame(value = TN_S/(TN_S+FP_S))
  TSS_S <- data.frame(value = SENSI_S[1] + SPECI_S[1] - 1)
  AUC_S <- data.frame(value = auc(obs_pred_S_cont$observed, obs_pred_S_cont$predicted))
  
  # stack the results of TSS and AUC for the five evaluation runs!
  TSS_S_5 <- rbind(TSS_S_5,TSS_S)
  AUC_S_5 <- rbind(AUC_S_5,AUC_S)
}

# Create data frame with the results of TSS evaluation (=Mean of 5 runs)
TSS_S <-  data.frame(value = mean(TSS_S_5$value), metric= "TSS", scale = "S", species = spec1)

# Create data frame with the results of AUC evaluation (=Mean of 5 runs)
AUC_S <- data.frame(value = mean(AUC_S_5$value), metric= "AUC", scale = "S", species = spec1)

# Evaluation with the Boyce-index an presence only evaluation metric (1 run)
boyce_S <- ecospat.boyce(prediction_continuous, S, nclass = 0, window.w="default", res=1000, PEplot = FALSE,method = 'spearman')
boyce_S <- boyce_S$cor
# Create data frame with the results of Boyce-index
boyce_S <- data.frame(value = boyce_S, metric= "BOYCE", scale = "S", species = spec1)


# Create data frame with the the results of TSS, AUC and Boyce-index
metric_S <- rbind(TSS_S,AUC_S,boyce_S)

### 3) Evaluation at spatial extent of the Racherin
# Extract binary predictions at extent of the G
pres_G <- data.frame(predicted= raster::extract(prediction_binary, G))
pres_G <- data.frame(observed = 1, predicted =pres_G$predicted)

# Extract continuous predictions at extent of the G
pres_G_cont <- data.frame(predicted= raster::extract(prediction_continuous, G))
pres_G_cont <- data.frame(observed = 1, predicted =pres_G_cont$predicted)

# Extract True positive and False negative binary predictions
TP_G <- sum(pres_G$predicted == 1)/nrow(pres_G)
FN_G <- sum(pres_G$predicted == 0)/nrow(pres_G)

# Prepare dataframes to store evaluation according to TSS and AUC
TSS_G_5 <- data.frame(value = as.numeric())
AUC_G_5 <- data.frame(value = as.numeric())
# Evaluation of TSS and AUC with five iterations
for (t in 1:5){
  # Sample random pseudo-absences within the layer created before
  absent <- data.frame(x=numeric(), y=numeric(), layer=numeric())
  # Sample as many pseudo-absences as species obserations used for model calibration
  for (i in 1:nrow(G)){
    # Create a 20 km buffer around each presence observation
    point <- G %>% slice(i) 
    coordinates(point) <- ~X_round + Y_round
    crs(point) <- crs(prediction_binary)
    point <- spTransform(point, crs(prediction_binary))
    point <- raster::buffer(point, width=20000)
    
    # Crop buffer to raster which excludes presences
    pot_ab <- crop(pot_absences, point)
    pot_ab <- mask(pot_ab, point)
    
    # Sample pseudo-absences where there are no presence observation (= !NA)
    absence <- data.frame(sampleRandom(pot_ab, 1, na.rm=T, cells=F, xy=T))
    absent <- rbind(absent,absence)
  }
  # Extract coordinates of sampled pseudo-absences
  absences<- absent[1:2]
  
  # Extract binary predictions of the pseudo-absences
  abs_G <- data.frame(predicted = raster::extract(prediction_binary, absences))
  abs_G <- data.frame(observed = 0, predicted =abs_G$predicted)
  
  # Extract binary predictions of the pseudo-absences
  abs_G_cont <- data.frame(predicted = raster::extract(prediction_continuous, absences))
  abs_G_cont <- data.frame(observed = 0, predicted =abs_G_cont$predicted)
  
  # Extract False positive and True negative binary predictions
  FP_G <- sum(abs_G$predicted == 1)/nrow(abs_G)
  TN_G <- sum(abs_G$predicted == 0)/nrow(abs_G) 
  
  # Create data frame with observed vs predicted values (Binary)
  obs_pred_G <- rbind(pres_G,abs_G)
  
  # Create data frame with observed vs predicted values (Continuous)
  obs_pred_G_cont <- rbind(pres_G_cont,abs_G_cont)
  
  
  # Calculate Sensitivity, Specificity, True skill statistic and Area under the curve at extent of the G
  SENSI_G <- data.frame(value = TP_G/(TP_G+FN_G))
  SPECI_G <- data.frame(value = TN_G/(TN_G+FP_G))
  TSS_G <- data.frame(value = SENSI_G[1] + SPECI_G[1] - 1)
  AUC_G <- data.frame(value = auc(obs_pred_G_cont$observed, obs_pred_G_cont$predicted))
  
  # stack the results of TSS and AUC for the five evaluation runs!
  TSS_G_5 <- rbind(TSS_G_5,TSS_G)
  AUC_G_5 <- rbind(AUC_G_5,AUC_G)
}

# Create data frame with the results of TSS evaluation (=Mean of 5 runs)
TSS_G <-  data.frame(value = mean(TSS_G_5$value), metric= "TSS", scale = "G", species = spec1)

# Create data frame with the results of AUC evaluation (=Mean of 5 runs)
AUC_G <- data.frame(value = mean(AUC_G_5$value), metric= "AUC", scale = "G", species = spec1)

# Evaluation with the Boyce-index an presence only evaluation metric (1 run)
boyce_G <- ecospat.boyce(prediction_continuous, G, nclass = 0, window.w="default", res=1000, PEplot = FALSE,method = 'spearman')
boyce_G <- boyce_G$cor
# Create data frame with the results of Boyce-index
boyce_G <- data.frame(value = boyce_G, metric= "BOYCE", scale = "G", species = spec1)


# Create data frame with the the results of TSS, AUC and Boyce-index
metric_G <- rbind(TSS_G,AUC_G,boyce_G)

### 4) Evaluation at spatial extent of the Hochschwab
# Extract binary predictions at extent of the H
pres_H <- data.frame(predicted= raster::extract(prediction_binary, H))
pres_H <- data.frame(observed = 1, predicted =pres_H$predicted)

# Extract continuous predictions at extent of the H
pres_H_cont <- data.frame(predicted= raster::extract(prediction_continuous, H))
pres_H_cont <- data.frame(observed = 1, predicted =pres_H_cont$predicted)

# Extract True positive and False negative binary predictions
TP_H <- sum(pres_H$predicted == 1)/nrow(pres_H)
FN_H <- sum(pres_H$predicted == 0)/nrow(pres_H)

# Prepare dataframes to store evaluation according to TSS and AUC
TSS_H_5 <- data.frame(value = as.numeric())
AUC_H_5 <- data.frame(value = as.numeric())
# Evaluation of TSS and AUC with five iterations
for (t in 1:5){
  # Sample random pseudo-absences within the layer created before
  absent <- data.frame(x=numeric(), y=numeric(), layer=numeric())
  # Sample as many pseudo-absences as species obserations used for model calibration
  for (i in 1:nrow(H)){
    # Create a 20 km buffer around each presence observation
    point <- H %>% slice(i) 
    coordinates(point) <- ~X_round + Y_round
    crs(point) <- crs(prediction_binary)
    point <- spTransform(point, crs(prediction_binary))
    point <- raster::buffer(point, width=20000)
    
    # Crop buffer to raster which excludes presences
    pot_ab <- crop(pot_absences, point)
    pot_ab <- mask(pot_ab, point)
    
    # Sample pseudo-absences where there are no presence observation (= !NA)
    absence <- data.frame(sampleRandom(pot_ab, 1, na.rm=T, cells=F, xy=T))
    absent <- rbind(absent,absence)
  }
  # Extract coordinates of sampled pseudo-absences
  absences<- absent[1:2]
  
  # Extract binary predictions of the pseudo-absences
  abs_H <- data.frame(predicted = raster::extract(prediction_binary, absences))
  abs_H <- data.frame(observed = 0, predicted =abs_H$predicted)
  
  # Extract binary predictions of the pseudo-absences
  abs_H_cont <- data.frame(predicted = raster::extract(prediction_continuous, absences))
  abs_H_cont <- data.frame(observed = 0, predicted =abs_H_cont$predicted)
  
  # Extract False positive and True negative binary predictions
  FP_H <- sum(abs_H$predicted == 1)/nrow(abs_H)
  TN_H <- sum(abs_H$predicted == 0)/nrow(abs_H) 
  
  # Create data frame with observed vs predicted values (Binary)
  obs_pred_H <- rbind(pres_H,abs_H)
  
  # Create data frame with observed vs predicted values (Continuous)
  obs_pred_H_cont <- rbind(pres_H_cont,abs_H_cont)
  
  
  # Calculate Sensitivity, Specificity, True skill statistic and Area under the curve at extent of the H
  SENSI_H <- data.frame(value = TP_H/(TP_H+FN_H))
  SPECI_H <- data.frame(value = TN_H/(TN_H+FP_H))
  TSS_H <- data.frame(value = SENSI_H[1] + SPECI_H[1] - 1)
  AUC_H <- data.frame(value = auc(obs_pred_H_cont$observed, obs_pred_H_cont$predicted))
  
  # stack the results of TSS and AUC for the five evaluation runs!
  TSS_H_5 <- rbind(TSS_H_5,TSS_H)
  AUC_H_5 <- rbind(AUC_H_5,AUC_H)
}

# Create data frame with the results of TSS evaluation (=Mean of 5 runs)
TSS_H <-  data.frame(value = mean(TSS_H_5$value), metric= "TSS", scale = "H", species = spec1)

# Create data frame with the results of AUC evaluation (=Mean of 5 runs)
AUC_H <- data.frame(value = mean(AUC_H_5$value), metric= "AUC", scale = "H", species = spec1)

# Evaluation with the Boyce-index an presence only evaluation metric (1 run)
boyce_H <- ecospat.boyce(prediction_continuous, H, nclass = 0, window.w="default", res=1000, PEplot = FALSE,method = 'spearman')
boyce_H <- boyce_H$cor
# Create data frame with the results of Boyce-index
boyce_H <- data.frame(value = boyce_H, metric= "BOYCE", scale = "H", species = spec1)


# Create data frame with the the results of TSS, AUC and Boyce-index
metric_H <- rbind(TSS_H,AUC_H,boyce_H)

# Create a data frame combining all evaluation results
metric_all <- rbind(metric_alps,metric_S,metric_G,metric_H)
# Combine evaluation metrics for all species
metric_final <- rbind(metric_final,metric_all)
}

# Write evaluation results as a csv.
write.csv(metric_final, "Evaluation_Disk20km_final.csv")