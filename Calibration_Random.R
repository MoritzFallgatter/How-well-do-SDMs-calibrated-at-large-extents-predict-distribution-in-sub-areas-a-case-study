####### Script for using biomod2 for modeling multiple species #########

# Calibrating SDMs with random PA-selection within the study area (= Alps + 20-kilometer buffer)

args = commandArgs(trailingOnly=TRUE)
myRespName <- args[1] 

# Define a set of predictor variables to be used in the modeling process
predi <- c("Bio4", "Bio10", "Bio18", "Bio19", "HLI", "TWI")

## Installing and loading necessary libraries
library(biomod2)  
library(terra)    
library(sp)       

#############################################
########## Loading Data ######################
#############################################

# Print message to indicate which species is being modeled
print(paste(myRespName, " modeling...", sep=""))

# Load species occurrence data from a CSV file for the specified species
# The CSV file should contain species occurrence data with columns including species coordinates
csvfile <- sprintf("../Input_data_all_alps_mountains/%s_all.csv", myRespName)
data_all <- read.csv(csvfile, sep=";", dec=",")  # Load species data
data <- subset(data_all, mountains == "Alps")   # Subset data for the Alps region

#############################################
######## Defining Biomod2 Input ##############
#############################################

# Modify species name for use in the model
myRespName <- gsub("-", ".", myRespName)

# Split the data into training (80%) and evaluation (20%) sets for model validation
index_param <- sample(1:nrow(data), floor(0.80 * nrow(data)))  # Randomly select 80% of the data
data_save <- data  # Create a copy of the data
data_save$param <- 0
data_save$param[index_param] <- 1  # Assign 1 to the training set
write.csv(data_save, paste0("data_model/", myRespName, ".csv"))  # Save the split data
data_param <- data[index_param,]  # Extract training data
data_eval <- data[!index_param,]  # Extract evaluation data

# Prepare the presence/absence data for the species
myResp <- as.numeric(data_param[,"occ"])  # Extract occurrence data for the species
myResp[] <- 1  # Set all presence values to 1
numOfPres <- length(myResp)  # Calculate the number of presence points

# Extract the coordinates of species occurrences (XY coordinates)
myRespCoord <- data_param[c("X_round", "Y_round")]

## Loading predictor variables 
infiles <- paste("../Predictors_alps_20km_buffer/", c("bio4_CHELSA__35a_ts.tif", "bio10_CHELSA__35a_ts.tif",
                                                      "bio18_CHELSA__35a_ts.tif", "bio19_CHELSA__35a_ts.tif", "HLI.tif", "Alps_TWI_100m.tif"), sep="")

myExpl <- rast(infiles)  # Load the raster files
names(myExpl) <- c("Bio4", "Bio10", "Bio18", "Bio19", "HLI", "TWI")  # Name the predictor variables

#############################################
######## 1. Formatting Data ##################
#############################################

# Format the data for use in BIOMOD2
myBiomodData <- BIOMOD_FormatingData(resp.var   = myResp,
                                     expl.var   = myExpl,
                                     resp.xy    = myRespCoord,
                                     resp.name  = myRespName,
                                     PA.nb.rep  = 5,  # Repetitions of PA selection
                                     PA.nb.absences = numOfPres,  # Number of background points = number of presence point
                                     PA.strategy = 'random',  # Strategy for generating background data
                                     na.rm       = TRUE)

#############################################
#### 2. Define Modeling Options ############
#############################################

# Define the set of modeling algorithms to be used
modelset <- c("RF", "GBM", "ANN")  # Random Forest, Generalized Boosted Models, and Artificial Neural Networks

# Set the modeling options
myBiomodOption <- BIOMOD_ModelingOptions()

#############################################
#### 3. Species Distribution Model (SDM) Calibration ##################
#############################################

#############################################
#### 3.1. Machine Learning ##################
#############################################

# Fit the models using BIOMOD2 with the specified options
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    models = modelset,
                                    bm.options = myBiomodOption,
                                    CV.nb.rep = 3,  # Number repetitions
                                    CV.perc = 0.80,  # Percentage of data used for training
                                    var.import = 3,  # Compute variable importance
                                    metric.eval = c("TSS", "ROC", "KAPPA"),  # Metrics used to evaluate model performance
                                    scale.models = TRUE, 
                                    CV.do.full.models = FALSE,
                                    modeling.id = myRespName)

# Retrieve model evaluation results
eval <- get_evaluations(myBiomodModelOut)
print(names(eval))  # Print evaluation metrics
print(eval)  # Print model evaluations

# Filter and select models with high performance (TSS > 0.7)
eval2 <- eval[complete.cases(eval$validation),]  # Remove models with missing evaluation data
eval2 <- subset(eval2, metric.eval == "TSS")  # Filter for TSS metric
models_good7 <- sum(eval2$validation > 0.7)  # Count models with TSS > 0.7
models_good6 <- sum(eval2$validation > 0.6)  # Count models with TSS > 0.6

# Get variable importance from the models
varimp <- get_variables_importance(myBiomodModelOut)

# Save evaluation and variable importance results to CSV files
write.table(eval, paste("stats/eval_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")
write.table(varimp, paste("stats/varimp_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")

# Define threshold for model selection based on evaluation scores
if(models_good7 > 4) {
  threshy <- 0.7
} else if(models_good6 > 4) {
  threshy <- 0.6
} else {
  threshy <- 0.5
}

# Ensemble modeling based on selected models
myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = 'all',
  em.by = 'all',
  em.algo = c("EMwmean"),  # Ensemble algorithm: weighted mean
  metric.eval = c('TSS'),
  metric.select.thresh = threshy,  # Use selected threshold
  metric.select = c("TSS"),
  var.import = 3) 

# Save ensemble model evaluation results to CSV
eval_EM <- get_evaluations(myBiomodEM)
write.table(eval_EM, paste("stats/ensemble_", myRespName, ".csv", sep=""), dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE, sep=";")

# Save ensemble model variable importance results
write.table(get_variables_importance(myBiomodEM), paste("stats/", myRespName, "_EM_var_importance.txt", sep=""), sep=",")

###########################
# 5. Making Projections ###
###########################

# Define the name for the current projection
projection_name = paste("current_", sep="")
print(paste("Modelling: ", projection_name, sep=""))

# Perform projection of the ensemble model to the current environmental data
projection <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                new.env = myExpl,
                                proj.name = projection_name,
                                models.chosen = "all",
                                metric.binary = "TSS",
                                on_0_1000 = TRUE,  # Scale projections to range 0-1000
                                compress = TRUE,
                                build.clamping.mask = FALSE)

# Generate ensemble forecasting based on the projections
myForecast <- BIOMOD_EnsembleForecasting(
  bm.proj = projection,
  bm.em = myBiomodEM,
  models.chosen = "all",
  metric.binary = "TSS",
  compress = TRUE)

# Print completion message
print(paste("Done", sep=" "))
