# How well do SDMs calibrated at large extents predict distribution in sub-areas: a case study
This repository contains data and R scripts used to explore how the accuracy of species distribution models (SDMs) changes when models calibrated at a large spatial extent (the entire European Alps) are used to make projections at smaller spatial scales (individual mountain ranges).
In detail, we explored the potential change of accuracy of SDMs calibrated at the extent of the whole European Alps when used for projections onto individual mountain ranges. We first related the occurrences of 16 alpine plant species across the European Alps to six topo-climatic predictors at a spatial resolution of 100 x 100 meters. We then projected the distribution of these species throughout the Alps, including three mountains with a particularly high number of species occurrences. Subsequently, we compared the accuracy of the model projections at the extent of the Alps to those at each of the individual mountains. Furthermore, we explored whether the differences in projection accuracy might be driven by differences in the distribution of species along environmental gradients between the calibration extent and the small areas, i.e., by differences in their realized niches.

The provided species occurrences and following R-scripts were used in this context:
Calibration_Disk20.R = Calibrating SDMs selecting pseudo-absences within a 20-kilometer buffer around the species occurrences
Calibration_Random.R = Calibrating SDMs with random PA-selection within the study area (= Alps + 20-kilometer buffer)
Data_Preparation.R = Prepare species occurrence data
Evaluation_Disk20km.R = Evaluating SDMs selecting pseudo-absences within a 20-kilometer buffer around the species occurrences
Evaluation_Random.R = Evaluating SDMs with random PA-selection within the study area (= Alps + 20-kilometer buffer)
OLS_regression.R = OLS regression Model: Drop in model accuracy ~ Sum of weighted D
Predictor_correlation.R = Intercorrelation of predictor variables
Predictors_similarity_Kol_smir_D.R = Calculating of Kolmogorov Smirnov D statistic between predictors at the extent of the Alps and each individual mountain
Species_occurrences.csv = Cleaned species occurrences used in this study
Speciesnames.csv = Names of the species used in this study
Sum_of_weighted_D.R = Kolmogorov Smirnov D statistic weighted by predictor importance and summed for each mountain (How dissimilar is each mountain to extent of Alps)
Wilcox_test.R = Signed Wilcox rank test between model accuracy at the extent of the Alps and each individual mountain for each evaluation metric (TSS, AUC and Boyce index), separately.
