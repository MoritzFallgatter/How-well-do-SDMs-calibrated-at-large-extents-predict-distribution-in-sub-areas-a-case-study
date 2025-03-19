# Correlation of Predictor Variables

# Clear the working environment
rm(list=ls())

# Load required libraries
library(rgdal)
library(raster)
library(tidyverse)

# Set working directory
setwd("C:/Users/mofal/Desktop/R-Publication/Data")

# Load species names 
# The file contains a column named 'species' with species names
species <- read.csv("Speciesnames.csv")
species <- species$species  # Extract species names into a vector
spec1 <- species[1]  # Select the first species for further analysis

# Load environmental predictor variables as a raster stack
predictors <- stack("Predictors/bio4_CHELSA__35a.ts.tif",
                    "Predictors/bio10_CHELSA__35a.ts.tif",
                    "Predictors/bio18_CHELSA__35a.ts.tif",
                    "Predictors/bio19_CHELSA__35a.ts.tif",
                    "Predictors/Alps_TWI_100m.tif",
                    "Predictors/HLI.tif")

# Rename raster layers for easier reference
names(predictors) <- c("Bio4", "Bio10", "Bio18", "Bio19", "TWI", "HLI")

# Create an empty data frame to store correlation results
result <- data.frame()

# Iterate over each pair of predictor variables to compute the correlation
# This nested loop calculates Pearson's correlation coefficient for each pair
for (i in 1:6){
  x <- values(predictors[[i]])  # Extract raster values for the i´th predictor
  
  for (t in i:6){
    y <- values(predictors[[t]])  # Extract raster values for the t´th predictor
    
    # Perform Pearson's correlation test between the two predictors
    correlation <- cor.test(x, y, method = 'pearson')
    
    # Round the correlation coefficient and p-value for readability
    cor <- round(correlation$estimate, 3)
    p <- round(correlation$p.value, 3)
    
    # Store the results in a data frame with the predictor variable names, correlation, and p-value
    res <- data.frame(x = names(predictors[[i]]),
                      y = names(predictors[[t]]),
                      cor = cor, p = p)
    
    # Stack the result to the data frame created before the loop
    result <- rbind(result, res)
  }
}

# Write the results to a CSV file 
write.csv(result, "Correlation_predictors.csv")