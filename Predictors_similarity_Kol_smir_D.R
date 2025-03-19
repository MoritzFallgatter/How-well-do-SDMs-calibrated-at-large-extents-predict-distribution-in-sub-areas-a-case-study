# Calculate the dissimilarity of predictors at the extent of the Alps to each of the individual mountains

# Clear the working environment
rm(list=ls())

# Load necessary libraries
library(rgdal)
library(raster)
library(tidyverse)

# Set working directory
setwd("C:/Users/mofal/Desktop/R-Publication/Data")

# Load species names from a CSV file
species <- read.csv("Speciesnames.csv")
species <- species$species # Extract the species column as a vector
spec1 <- species[1] # Select the first species for initialization

# Load environmental predictor variables as raster stack
predictors <- stack("Predictors/bio4_CHELSA__35a.ts.tif",
                    "Predictors/bio10_CHELSA__35a.ts.tif",
                    "Predictors/bio18_CHELSA__35a.ts.tif",
                    "Predictors/bio19_CHELSA__35a.ts.tif",
                    "Predictors/Alps_TWI_100m.tif",
                    "Predictors/HLI.tif")

# Rename the layers of the predictor stack for better readability
names(predictors) <- c("Bio4", "Bio10", "Bio18", "Bio19", "TWI", "HLI")

# Initialize an empty data frame to store Kolmogorov-Smirnov test results (= "D")
result <- data.frame(alpen_SK = as.numeric(), alpen_GK = as.numeric(), alpen_HS = as.numeric(), 
                     predictor = as.character(), species = as.character(), Pa_approach = as.character())

# Define the model approaches (random and disk)
model <- c("random", "disk")

# Loop through each species for analysis
for(spec1 in species){
  
  # Modify species names to match the data files (replacing "-" with ".")
  spec <- gsub("-",".",spec1, fixed=TRUE)
  
  # Loop through both model approaches
  for(y in 1:2){
    # Define the file path for species observations within the extent of the Alps
    input_file_alps <- paste("moritz_johannes_20240416/", model[y], "/data_model/", spec, ".csv", sep="")
    
    # Load species observation data and convert coordinate values to numeric
    input_data_alps <- read.csv(input_file_alps, sep=",")
    input_data_alps$X_round <- as.numeric(gsub(",", ".", input_data_alps$X_round, fixed=TRUE))
    input_data_alps$Y_round <- as.numeric(gsub(",", ".", input_data_alps$Y_round, fixed=TRUE))
    
    # Define the file path for species observations at the individual mountains
    input_file_mountains <- paste("Input_data_all_alps/", spec1, "_all", ".csv", sep="")
    
    # Load species occurrences at the individual mountains and convert coordinate values to numeric
    input_file_mountains <- read.csv(input_file_mountains, sep=";")
    input_file_mountains$X_round <- as.numeric(gsub(",", ".", input_file_mountains$X_round, fixed=TRUE))
    input_file_mountains$Y_round <- as.numeric(gsub(",", ".", input_file_mountains$Y_round, fixed=TRUE))
    
    # Filter species observations which were used for model calibration at the extent of the Alps
    alps <- input_data_alps %>% filter(param == 1) %>%
      dplyr::select(X_round, Y_round) 
    
    # Filter species observations for specific mountain ranges (Schrankogel, Racherin, Hochschwab)
    S <- input_file_mountains %>% filter(mountains == "S") %>%
      dplyr::select(X_round, Y_round) 
    G <- input_file_mountains %>% filter(mountains == "G") %>%
      dplyr::select(X_round, Y_round) 
    H <- input_file_mountains %>% filter(mountains == "H") %>%
      dplyr::select(X_round, Y_round) 
    
    # Loop through each environmental predictor variable
    for(i in 1:6){
      
      # Extract the environmental conditions from each raster layer at each extent
      alpen <- raster::extract(predictors[[i]], alps)
      SK    <- raster::extract(predictors[[i]], S)
      GK    <- raster::extract(predictors[[i]], G)
      HS    <- raster::extract(predictors[[i]], H)
      
      # Perform the Kolmogorov-Smirnov test between the environmental conditions at the extent of the Alps and each indivudial mountain
      # The D statistic measures the dissimilarity between distributions, with smaller values indicating more similarity
      alpen_SK <- ks.test(alpen, SK)
      alpen_GK <- ks.test(alpen, GK)
      alpen_HS <- ks.test(alpen, HS)
      
      # Store the results of the KS test for each predictor variable and model approach
      res <- data.frame(alpen_SK[1], alpen_GK[1], alpen_HS[1],
                        names(predictors[[i]]), spec, model[y])
      
      # Rename columns for clarity
      colnames(res) <- c("Schrankogel", "Racherin", "Hochschwab", "Predictor", "Species", "Pa_approach")
      
      # Combine the results for all species and predictor variables
      result <- rbind(result, res) 
    }
  }
}

# Save the results to a CSV file for further analysis
write.csv(result, "Predictors_similarity_separated_by_pa.csv")
