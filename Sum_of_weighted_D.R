# Calculation of the Sum of weighted dissimilarity D 

# Clear the working environment
rm(list=ls())

# Load necessary libraries
library(tidyverse)
library(dplyr)
library(ggpubr)

# Set the working directory
setwd("C:/Users/mofal/Desktop/R-Publication/Data")

# Load the species names
species <- read.csv("Speciesnames.csv")
species <- species$species # Extract the species column as a vector
spec1 <- species[1] # Select the first species for initialization

# Define the model approaches and predictor variables
pa_approach <- c("random", "disk")
pred <- c("Bio4", "Bio10", "Bio18", "Bio19", "TWI", "HLI")

# Load Data
# Load the variable importance data
var_imp <- read.csv("Variable_importance_EM.csv")
var_imp <- var_imp[2:5] # Select the relevant columns
# Load the Kolmogorov-Smirnov test results for predictor similarity by PA approach
kol_smir <- read.csv2("Predictors_similarity_sepetated_by_pa.csv", sep = ",") 
kol_smir <- kol_smir[2:6] # Select the relevant columns

# List of mountain regions to process
mountain <- c("Schrankogel", "Hochschwab", "Racherin")

# Initialize empty data frames to store the results for each model approach
result_random <- data.frame()
result_disk20 <- data.frame()
result_all <- data.frame()

# Loop through each species to calculate the Sum of weighted D
for (spec1 in species) {
  
  # Modify species names to fit data file format (replace "-" with ".")
  spec <- gsub("-", ".", spec1, fixed = TRUE)
  
  # Filter variable importance and KS similarity results for the current species
  varimp <- filter(var_imp, Species == spec)
  kolsmir <- filter(kol_smir, Species == spec)
  # Flip the data for furhter processing
  kolsmir <- pivot_longer(kolsmir, 1:3, values_to = "D", names_to = "Mountain")
  
  # Multiply D statistic of each variable by the variable importance for each mountain
  # Merge variable importance and KS similarity data for the species
  data <- left_join(varimp, kolsmir, by = c("Species", "Predictor"))
  data$D <- as.numeric(data$D) # Convert the D statistic to numeric
  data$Weighted_D <- data$Varimp * data$D # Calculate the weighted D (importance * dissimilarity)
  
  # Separate the data by PA approach (random vs disk)
  random <- filter(data, Pa_approach == "random")
  disk20 <- filter(data, Pa_approach == "disk")
  
  # Sum up the weighted D values per mountain region and species for each PA approach
  # Look through mountains
  for (i in 1:3) {
    mount_random <- filter(random, Mountain == mountain[i])
    mount_disk <- filter(disk20, Mountain == mountain[i])
    
    # Calculate the mean weighted D for the random model approach
    meanDweighted_r <- sum(mount_random$Weighted_D)
    # Calculate the mean weighted D for the disk model approach
    meanDweighted_d <- sum(mount_disk$Weighted_D)
    
    # Store the results for the random approach
    result_rand <- data.frame(Species = spec1, Pa_approach = "random", Mountain = mountain[i], Weighted_D = meanDweighted_r)
    result_random <- rbind(result_random, result_rand)
    
    # Store the results for the disk approach
    result_disk <- data.frame(Species = spec1, Pa_approach = "disk", Mountain = mountain[i], Weighted_D = meanDweighted_d)
    result_disk20 <- rbind(result_disk20, result_disk)
  }
}

# Combine the results from both PA approaches into one final data frame
result_all <- rbind(result_random, result_disk20)

# Save the final results
write.csv(result_all, "Weighted_D.csv")
