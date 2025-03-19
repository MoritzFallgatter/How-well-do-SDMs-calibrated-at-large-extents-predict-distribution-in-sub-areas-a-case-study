######################################################################################
####### Script for using biomod2 to model multiple species #####
######################################################################################

# Script for SDMs selecting pseudo-absences within a 20-kilometer buffer around the species occurrences

args = commandArgs(trailingOnly=TRUE)
myRespName <- args[1]

# Define predictor variables for the model
predi <- c("Bio4","Bio10","Bio18","Bio19","HLI","TWI")

# Load the necessary libraries
library(biomod2)
library(terra)
library(sp)

#############################################
########## Loading the Data #################
#############################################

# Loading species occurrence and absence data
# The model will be run for a specific species based on the provided argument (myRespName)

print(paste(myRespName, " modeling...", sep=""))

# Loading the data for the specified species
csvfile <- sprintf("../Input_data_all_alps_mountains/%s_all.csv", myRespName)
data_all <- read.csv(csvfile, sep=";", dec=",")
# Subset data to only include presence observations at extent of the Alps
data <- subset(data_all, mountains == "Alps")

#############################################
######## Defining Biomod Input ##############
#############################################
# Change species name for use in file paths
myRespName <- gsub("-", ".", myRespName)

# Randomly sample data for model training and testing (80% training, 20% testing)
index_param <- sample(1:nrow(data), floor(0.80 * nrow(data)))
data_save <- data
data_save$param <- 0
data_save$param[index_param] <- 1

# Save the processed data with param column indicating training and testing sets
write.csv(data_save, paste0("data_model/", myRespName, ".csv"))

# Separate the data into training and evaluation sets
data_param <- data[index_param, ]
data_eval <- data[!index_param, ]

# Defining presence data for species occurrence
myResp <- as.numeric(data_param[,"occ"])
myResp[] <- 1
numOfPres <- length(myResp)  # Count of presence occurrences

# XY coordinates for species occurrences (2 columns: X and Y)
myRespCoord <- data_param[c("X_round", "Y_round")]

# Load predictor variables
infiles <- paste("../Predictors_alps_20km_buffer/", c("bio4_CHELSA__35a_ts.tif", "bio10_CHELSA__35a_ts.tif", "bio18_CHELSA__35a_ts.tif", "bio19_CHELSA__35a_ts.tif", "HLI.tif", "Alps_TWI_100m.tif"), sep="")
myExpl <- rast(infiles)
names(myExpl) <- c("Bio4", "Bio10", "Bio18", "Bio19", "HLI", "TWI")

# Formatting the data for Biomod2
myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,
  expl.var = myExpl,
  #eval.resp.var = myResp_eval,
  resp.xy = myRespCoord,
  #eval.resp.xy = myRespCoord_eval,
  resp.name = myRespName,
  PA.nb.rep = 5,  # Number of pseudo-absence repetitions
  PA.nb.absences = numOfPres,  # Number of pseudo-absences
  PA.strategy = 'disk',  # Strategy for pseudo-absence generation
  PA.dist.min = 0,  # Minimum distance for pseudo-absence selection
  PA.dist.max = 20000,  # Maximum distance for pseudo-absence selection
  na.rm = TRUE  # Remove NA values from the data
)

#############################################
#### 2. Define Modeling Options #############
#############################################

# Set the modeling algorithms to be used
modelset <- c("RF", "GBM", "ANN")

# Use default settings for modeling options
myBiomodOption <- BIOMOD_ModelingOptions()

#############################################
#### 3. Model Calibration ###################
#############################################

# Machine learning: Running the models
myBiomodModelOut <- BIOMOD_Modeling(
  bm.format = myBiomodData, 
  models = modelset,
  bm.options = myBiomodOption, 
  CV.nb.rep = 3,  # Number of repetitions
  CV.perc = 0.80,  # Proportion of data used for training (80%)
  var.import = 3,  # Compute variable importance
  metric.eval = c("TSS", "ROC", "KAPPA"),  # Evaluation metrics
  scale.models = TRUE,  
  CV.do.full.models = FALSE,  
  modeling.id = myRespName  
)

# Get evaluation results for all models
eval <- get_evaluations(myBiomodModelOut)
print(names(eval))
print(eval)  

# Filter and process evaluations for TSS metric
eval2 <- eval[complete.cases(eval$validation),]
eval2 <- subset(eval2, metric.eval == "TSS")

# Count the number of models with good performance (TSS > 0.7 or 0.6)
models_good7 <- sum(eval2$validation > 0.7)
models_good6 <- sum(eval2$validation > 0.6)

# Get variable importance for the models
varimp <- get_variables_importance(myBiomodModelOut)

# Save evaluation results and variable importance to files
write.table(eval, paste("stats/eval_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")
write.table(varimp, paste("stats/varimp_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")

# Define threshold based on model performance
if (models_good7 > 4) {
  threshy <- 0.7
} else if (models_good6 > 4) {
  threshy <- 0.6
} else {
  threshy <- 0.5
}

# Run ensemble modeling based on the selected threshold
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  em.algo = c("EMwmean"),
  metric.eval = c('TSS'),
  metric.select.thresh = threshy,
  metric.select = c("TSS"),
  var.import = 3
)

# Get evaluation results for ensemble model
eval_EM <- get_evaluations(myBiomodEM)

# Save ensemble model evaluation results
write.table(eval_EM, paste("stats/ensemble_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")

# Save variable importance for ensemble model
write.table(get_variables_importance(myBiomodEM), paste("stats/", myRespName, "_EM_var_importance.txt", sep=""), sep=",")

###########################
# 5. Making Projections ###
###########################

# Define projection name and run the projection
projection_name = paste("current_", sep="")
print(paste("Modeling: ", projection_name, sep=""))

projection <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut, 
  new.env = myExpl, 
  proj.name = projection_name, 
  models.chosen = "all", 
  metric.binary = "TSS", 
  on_0_1000 = TRUE, 
  compress = TRUE,  
  build.clamping.mask = FALSE 
)

# Run ensemble forecasting
myForecast <- BIOMOD_EnsembleForecasting(
  bm.proj = projection,
  bm.em = myBiomodEM,
  models.chosen = "all",
  metric.binary = "TSS",
  compress = TRUE
)

print(paste("Done", sep=" "))
