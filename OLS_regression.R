# OLS regression Model: Weighted Analysis for Model Accuracy Drop

# Clear working environment
rm(list=ls())

# Load necessary libraries
library(tidyverse) 
library(performance)
library(stats)  
library(ggpubr)

# Set the working directory where the data files are stored
setwd("C:/Users/mofal/Desktop/R-Publication/Data")

# Load the Weighted D statistics from a CSV file
Weighted_D <- read.csv("Weighted_D.csv")

# Define vectors to set up loops for various combinations of models and evaluation metrics
metric <- c("TSS", "AUC", "BOYCE")  # Evaluation metrics
random_disk20 <- c("Random_final", "Disk20km_final")  # Pseudo-absence selection approaches
modelname <- c("Random", "Disk20")  # Model names for reporting
modelname2 <- c("random", "disk")  # Alternative model names for filtering
predictors <- c("Bio4", "Bio10", "Bio18", "Bio19", "TWI", "HLI")  # Predictor variables

# Initialize empty data frames to store results
result <- data.frame(model = as.character(), metric = as.character(), predictor = as.character(), 
                     intercept = as.numeric(), estimate_p = as.numeric(), p_value = as.numeric(), 
                     R2_multiple = as.numeric(), R2_adjusted = as.numeric())

result_all <- data.frame(model = as.character(), metric = as.character(), predictor = as.character(), 
                         intercept = as.numeric(), estimate_p = as.numeric(), p_value = as.numeric(), 
                         R2_multiple = as.numeric(), R2_adjusted = as.numeric())

df_final <- data.frame()  # Data frame to store final combined results for model evaluation

# Loop through pseudo-absence selection approaches (Random and Disk20)
for(x in 1:2){
  
  # Select the current pseudo-absence approach and corresponding model name
  mo <- random_disk20[x]
  model_name <- modelname[x]
  pa_approach <- modelname2[x]
  
  # Load the model evaluation results (TSS, AUC, Boyce-index) for the selected pseudo-absence approach
  model <- read.csv2(paste("Evaluation_", mo, ".csv", sep=""), sep=",")
  
  # Filter the weighted D statistics for the selected PA approach
  Weighted_D_use <- filter(Weighted_D, Pa_approach == pa_approach)
  Weighted_D_use <- Weighted_D_use[2:5]  # Keep relevant columns
  
  # Loop through the evaluation metrics (TSS, AUC, Boyce) for the model evaluation
  for(i in 1:3){
    metric_use <- metric[i]  # Select the current evaluation metric
    
    # Filter the model evaluation results for the selected metric
    model_use <- model %>% dplyr::filter(metric == metric_use)
    model_use$value <- as.numeric(model_use$value)  # Convert values to numeric
    
    # For AUC, convert the values to the [-1, 1] range
    if(metric_use == "AUC"){
      model_use$value <- model_use$value * 2 - 1
    }
    
    # Flip the data frame for further procession
    model_use <- pivot_wider(model_use, names_from = scale)
    
    # Calculate the drop in model accuracy between the extent of the Alps and each individual mountain (Schrankogel, Racherin, Hochschwab)
    model_use <- data.frame(species = unique(model_use$species), 
                            Schrankogel = na.omit(model_use$ALPS) - na.omit(model_use$S),
                            Racherin = na.omit(model_use$ALPS) - na.omit(model_use$G), 
                            Hochschwab =  na.omit(model_use$ALPS) - na.omit(model_use$H))
    
    # Flip the data to long format for further analysis
    model_use <- pivot_longer(model_use, 2:4, names_to = "Mountain", values_to = "Drop")
    model_use$Metric <- metric_use  # Add the current evaluation metric to the data
    
    # Combine the model evaluation results with the weighted D statistics
    model_use2 <- left_join(Weighted_D_use, model_use)
    
    # Append the combined results to the final data frame
    df_final <- rbind(df_final, model_use2)
    model_use2$Species <- factor(model_use2$Species)  # Convert Species to a factor
    
    # Rename columns for clarity
    colnames(model_use2) <- c("Species", "Pa_approach", "Mountain", "Weighted_D", "Drop", "Metric")
    
    # Conduct linear regression: Evaluate the relationship between 'Drop' and 'Weighted_D'
    lm1 <- lm(Drop ~ Weighted_D , data = model_use2)
    summary <- summary(lm1)  # Get the summary of the linear model
    
    # Extract the estimates, p-values and R from the summary
    estimate <- data.frame(summary[4])
    r2 <- data.frame(summary[8])  # R-squared value
    r2adj <- data.frame(summary[9])  # Adjusted R-squared value
    
    # Prepare a result data frame with key model statistics and adjust p-values for multiple testing
    result <- data.frame(model = model_name, metric = metric_use, intercept = round(estimate[1,1], digits = 3),
                         std = round(estimate[1,2], digits = 3), p_value = p.adjust(estimate[2,4], method = "bonferroni"),
                         R2_multiple = round(r2[1,1], digits = 3), R2_adjusted = round(r2adj[1,1], digits = 3))
    
    # Reset temporary variables for the next loop iteration
    summary <- 0
    p_adj <- 0
    estimate <- 0
    
    # Append the result to the combined results data frame
    result_all <- rbind(result_all, result)
  }
}

# Save the final regression results to a CSV file
write.csv(result_all, "OLS_regression.csv")
