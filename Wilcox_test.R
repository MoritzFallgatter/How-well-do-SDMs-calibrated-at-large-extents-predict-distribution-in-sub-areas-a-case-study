# Wilcox-Test to investigating the differences in model accuracy between the spatial Extent of the Alps the the individual mountains

# Clear the environment
rm(list=ls())

# Load necessary libraries
library(ggpubr)
library(ggplot2)
library(tidyverse)
library(exactRankTests)

# Set working directory (adjust as needed)
setwd("C:/Users/mofal/Desktop/R-Publication/Data")

### Load model evaluation statistics
random <- read.csv("Evaluation_random_final.csv",sep=",")
random$model <- "Random"
disk20 <- read.csv("Evaluation_Disk20km_final.csv",sep=",")
disk20$model <- "Disk20"
data <- rbind(random, disk20)

# Build new dataframe with all the evaluation metrics included
tss <- filter(data, metric=="TSS")
auc <- filter(data, metric=="AUC")
# AUC adjusted to range from -1 to 1
auc$value <- auc$value*2-1
boyce <- filter(data, metric=="BOYCE")
data <- rbind(tss,auc,boyce)

# Convert evaluation metric values to numeric
# Replacing commas with periods for correct number formatting
data$value <- gsub(",",".",data$value, fixed =TRUE)
data$value <- as.numeric(data$value)

# Rename evaluation metrics for clarity
data$metric[data$metric == "BOYCE"] <-"Boyce index"
data$metric[data$metric == "AUC"] <- "AUC adjusted"

# Convert metrics to factors with specific order
data$metric = factor(data$metric, levels=c('TSS','AUC adjusted','Boyce index'))

# Standardize naming 
colnames(data)[colnames(data) == "scale"] <- "Extent"
colnames(data)[colnames(data) == "metric"] <- "Metric"
data$Extent[data$Extent == "ALPS"] <- "Alps"
data$Extent[data$Extent == "G"] <- "Racherin"
data$Extent[data$Extent == "S"] <- "Schrankogel"
data$Extent[data$Extent == "H"] <- "Hochschwab"

# Slit data by pseudo absence approach 
disk <- data %>% filter(model == "Disk20")
rand <- data %>% filter(model == "Random")

# Convert categorical variables to factors
rand$Extent <- factor(rand$Extent)
rand$Metric <- factor(rand$Metric)
rand$species <- factor(rand$species)
rand$model <- factor(rand$model)
rand$Extent <- factor(rand$Extent, levels = c("Alps", "Hochschwab", "Racherin", "Schrankogel"))
rand<-rand[2:5] # Select relevant columns

# Initialize results dataframe for Wilcoxon test
results_random <- data.frame()
met <- c("TSS","AUC adjusted","Boyce index")
mount <- c("Schrankogel", "Racherin", "Hochschwab")

# Perform Wilcoxon test for each metric and mountain
for(i in 1:3){
  mountain2 <- mount[i]
  data <- filter (rand, Extent %in% c("Alps", mountain2))
  for(x in 1:3){
    metric2 <- met[x]
    data2 <- filter(data, Metric == metric2)
    res <- exactRankTests::wilcox.exact(value ~ Extent, data = data2, paired = TRUE, exact=T)
    result <- data.frame(Extent = mountain2, Metric= metric2, Pvalue =as.numeric(res[3]), V = mean(data2$value))
    res<-0
    results_random <- rbind(results_random, result)
  }}

# Correct for multiple testing using Bonferroni adjustment
Schrankogel <- filter(results_random, Extent == "Schrankogel")
Schrankogel$adjusted_Pvalue <- p.adjust(Schrankogel$Pvalue, method = "bonferroni")
Racherin <- filter(results_random, Extent == "Racherin")
Racherin$adjusted_Pvalue <- p.adjust(Racherin$Pvalue, method = "bonferroni")
Hochschwab <- filter(results_random, Extent == "Hochschwab")
Hochschwab$adjusted_Pvalue <- p.adjust(Hochschwab$Pvalue, method = "bonferroni")

# Merge data back together
results_random <- rbind(Schrankogel, Racherin, Hochschwab)

# Assign significance levels based on adjusted p-values
for(y in 1:9){ 
  if (results_random$adjusted_Pvalue[y] <= 0.001){
    results_random$p[y]<- "***"
  }else if(results_random$adjusted_Pvalue[y] <= 0.01){
    results_random$p[y]<- "**"
  }else if(results_random$adjusted_Pvalue[y] <= 0.05){
    results_random$p[y]<- "*"
  }else{results_random$p[y]<- "n.s."}}

# Standardize factor levels
results_random$Extent <- factor(results_random$Extent, levels = c("Alps", "Hochschwab", "Racherin", "Schrankogel"))
results_random$Metric = factor(results_random$Metric, levels=c('TSS','AUC adjusted','Boyce index'))


##################  
# Standardize naming
disk$Extent <- factor(disk$Extent)
disk$Metric <- factor(disk$Metric)
disk$species <- factor(disk$species)
disk$model <- factor(disk$model)
disk$Extent <- factor(disk$Extent, levels = c("Alps", "Hochschwab", "Racherin", "Schrankogel"))
disk<-disk[2:5]

# Initialize results dataframe for Wilcoxon test
results_disk20 <- data.frame()
met <- c("TSS","AUC adjusted","Boyce index")
mount <- c("Schrankogel", "Racherin", "Hochschwab")

# Perform Wilcoxon test for each metric and mountain
for(i in 1:3){
  mountain2 <- mount[i]
  data <- filter (disk, Extent %in% c("Alps", mountain2))
  for(x in 1:3){
    metric2 <- met[x]
    data2 <- filter(data, Metric == metric2)
    res <- exactRankTests::wilcox.exact(value ~ Extent, data = data2, paired = TRUE, exact=T)
    result <- data.frame(Extent = mountain2, Metric= metric2, Pvalue =as.numeric(res[3]), V = mean(data2$value))
    res<-0
    results_disk20 <- rbind(results_disk20, result)
  }}

# Correct for multiple testing using Bonferroni adjustment
Schrankogel <- filter(results_disk20, Extent == "Schrankogel")
Schrankogel$adjusted_Pvalue <- p.adjust(Schrankogel$Pvalue, method = "bonferroni")
Racherin <- filter(results_disk20, Extent == "Racherin")
Racherin$adjusted_Pvalue <- p.adjust(Racherin$Pvalue, method = "bonferroni")
Hochschwab <- filter(results_disk20, Extent == "Hochschwab")
Hochschwab$adjusted_Pvalue <- p.adjust(Hochschwab$Pvalue, method = "bonferroni")

# Merge data back together
results_disk20 <- rbind(Schrankogel, Racherin, Hochschwab)

# Assign significance levels based on adjusted p-values
for(y in 1:9){ 
  if (results_disk20$adjusted_Pvalue[y] <= 0.001){
    results_disk20$p[y]<- "***"
  }else if(results_disk20$adjusted_Pvalue[y] <= 0.01){
    results_disk20$p[y]<- "**"
  }else if(results_disk20$adjusted_Pvalue[y] <= 0.05){
    results_disk20$p[y]<- "*"
  }else{results_disk20$p[y]<- "n.s."}}

# Standardize factor levels
results_disk20$Extent <- factor(results_disk20$Extent, levels = c("Alps", "Hochschwab", "Racherin", "Schrankogel"))
results_disk20$Metric = factor(results_disk20$Metric, levels=c('TSS','AUC adjusted','Boyce index'))

# Define theme for plots
mytheme <- theme_classic()+
  theme(axis.text.x = element_text(size=18, angle = 15, hjust=0.5, vjust= 0.5), axis.text.y = element_text(size=18),title = element_text(size=18), axis.title = element_text(size=18), strip.text = element_text(size = 18),panel.spacing = unit(2, "lines"))

# Create and save plots
ggrand <- ggplot(rand, aes(x=Extent, y=value, fill =Extent))+
  facet_wrap(~Metric)+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("white","#1b9e77", "#d95f02", "#7570b3"))+
  geom_jitter(width=0.2,shape = 21,size=3)+
  geom_text(data = results_random, aes(x=Extent, y=1.1, label = p), size=7)+
  mytheme+  
  theme(legend.position = "none")+
  labs(title="", x= "",y="Model accuracy", legend="")+ylim(-1,1.1)

windows()
ggrand

# Save plot
ggsave("Plots_final/Evaluation_Random.png", width = 375, height = 200 , units = "mm")

# Create and save plots
ggdisk <- ggplot(disk, aes(x=Extent, y=value, fill =Extent))+
  facet_wrap(~Metric)+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("white","#1b9e77", "#d95f02", "#7570b3"))+
  geom_jitter(width=0.2,shape = 21,size=3)+
  geom_text(data = results_disk20, aes(x=Extent, y=1.1, label = p), size=7)+
  mytheme+
  theme(legend.position = "none")+
  labs(title="", x= "",y="Model accuracy", legend="")+ ylim(-1,1.1)

windows()
ggdisk

# Save plot
ggsave("Plots_final/Evaluation_Disk20.png", width = 375, height = 200 , units = "mm")

