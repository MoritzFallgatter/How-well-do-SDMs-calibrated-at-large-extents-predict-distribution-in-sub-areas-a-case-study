#### Prepare species occurrence data

# Load required libraries
library(raster) 
library(rgdal)  
library(sp)      
library(CoordinateCleaner)  

# Set the working directory
setwd("C:/Users/mofal/Desktop/Uni/master")

# Define the resolution of the raster grid (grid cell size for aggregating points)
resolution <- 100

# Define the spatial extent of the study area using a raster file
clim <- raster("alps_buffer_20km.tiff")  # Load the climate data raster
ext <- extent(clim)  # Get the extent (bounding box) of the raster

# Read the species occurrence data from the Austrian Vegetation Database
AT_Database <- read.csv("Data_ForMoritz.csv", sep = ";")

# Sort the species names alphabetically and create a unique list of species
species_w <- data.frame(species = sort(unique(AT_Database$Taxon_name.MC)))
dim(species_w)  # Display dimensions to verify: 233 species recorded

# Extract species list from available binary files for GBIF
species_g <- list.files("C:/Users/mofal/Desktop/master/SpeciesData", pattern =".bin$")
species_g <- gsub(".bin", "", species_g)  # Remove ".bin" extension from filenames
length(species_g)  # Display number of species recorded in GBIF: 249 species

# Identify species present only in one dataset (AT_Database or GBIF)
species_onlyW <- species_w[!species_w[,1] %in% species_g, 1]  # Species only in AT_Database
species_onlyG <- species_g[!species_g %in% species_w[,1]]  # Species only in GBIF

# Set the species list to process for further analysis (using GBIF species here)
species <- species_g

# Initialize an empty data frame for storing statistics
all_stats <- NULL 
# Loop through each species to process their occurrence data
for(spec in species) {
  print(spec)  # Print the current species name
  name_data <- gsub(" ", "-", spec)  # Replace spaces in species name with dashes
  name_w <- spec  # Species name for AT_Database
  name_g <- spec  # Species name for GBIF
  
  # Special case handling for taxonomic correction
  if(name_w == "Anthyllis vulneraria subsp. alpicola"){
    name_w <- "Anthyllis vulneraria"
  }
  
  ## Processing AT_Database's data
  data_w <- subset(AT_Database, Taxon_name.MC == name_w)  # Subset data for the current species
  if(nrow(data_w) > 0) {  # If there is data for the species in the AT_Database
    # Convert coordinates to spatial points and transform projection to EPSG:3035
    coordinates(data_w) = ~longitude + latitude
    proj4string(data_w) <- crs("EPSG:4326")  # Set WGS 84 coordinate reference system
    data_laea <- spTransform(data_w, crs("EPSG:3035"))  # Project to Lambert Azimuthal Equal Area
    data_laea <- data.frame(data_laea)  # Convert back to a data frame
    data_laea$long_round <- floor(data_laea$longitude / resolution) * resolution + resolution / 2
    data_laea$lat_round <- floor(data_laea$latitude / resolution) * resolution + resolution / 2
    data_laea$occ <- 1  # Set occurrence indicator
    # Aggregate occurrences by grid cells
    res_w <- aggregate(data_laea$occ, by = list(X_round = data_laea$long_round, Y_round = data_laea$lat_round), FUN = sum)
    names(res_w)[3] <- name_data  # Rename the third column to the species name
    write.csv2(res_w, paste("occurrences_cleaned/", name_data, "_W.csv", sep = ""), row.names = FALSE)  # Save to CSV
  } else {
    res_w <- data_w  # If no data, return original data (empty)
  }
  
  ## Processing GBIF data
  filename <- paste("SpeciesData/", spec, ".bin", sep = "")  # Path to GBIF binary data file
  load(filename)  # Load the GBIF data
  data_g <- data  # Store GBIF data
  
  # Filter GBIF data by year and coordinate uncertainty
  data_time <- subset(data_g, year >= 1980 & year <= 2022)  # Filter for recent years
  data_certain <- subset(data_time, coordinateUncertaintyInMeters <= resolution)  # Filter by uncertainty threshold
  data_certain$decimalLongitude <- as.numeric(data_certain$decimalLongitude)  # Convert longitude to numeric
  data_certain$decimalLatitude <- as.numeric(data_certain$decimalLatitude)  # Convert latitude to numeric
  
  if(nrow(data_certain) > 0) {  # If there is filtered data for the species in GBIF
    # Clean coordinates using predefined tests (e.g., checking for valid locations)
    flags <- clean_coordinates(x = data_certain, lon = "decimalLongitude", lat = "decimalLatitude",
                               species = "species", tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros"))
    data_cl <- data_certain[flags$.summary, ]  # Keep only valid records
    coordinates(data_cl) = ~decimalLongitude + decimalLatitude  # Convert to spatial points
    proj4string(data_cl) <- crs("EPSG:4326")  # Set WGS 84 CRS
    data_cl_laea <- spTransform(data_cl, crs("EPSG:3035"))  # Transform to Lambert Azimuthal Equal Area projection
    data_cl_laea <- data.frame(data_cl_laea)  # Convert back to data frame
    data_cl_laea$long_round <- floor(data_cl_laea$decimalLongitude / resolution) * resolution + resolution / 2
    data_cl_laea$lat_round <- floor(data_cl_laea$decimalLatitude / resolution) * resolution + resolution / 2
    data_cl_laea$occ <- 1  # Set occurrence indicator
    res_g1 <- aggregate(data_cl_laea$occ, by = list(X_round = data_cl_laea$long_round, Y_round = data_cl_laea$lat_round), FUN = sum)
    res_g <- subset(res_g1, X_round > ext@xmin & X_round < ext@xmax & Y_round > ext@ymin & Y_round < ext@ymax)  # Filter by study area extent
    names(res_g)[3] <- name_data  # Rename the third column to species name
    write.csv2(res_g, paste("occurrences_cleaned/", name_data, "_G.csv", sep = ""), row.names = FALSE)  # Save to CSV
  } else {
    res_g <- data_certain  # If no data, return original data (empty)
  }
  
  # Combine occurrences from both datasets (AT_Database and GBIF)
  data_all <- rbind(res_w, res_g)
  names(data_all)[3] <- "occ"  # Rename column to "occ" for occurrence count
  res_all <- aggregate(data_all$occ, by = list(X_round = data_all$X_round, Y_round = data_all$Y_round), FUN = sum)
  names(res_all)[3] <- name_data  # Rename the third column to species name
  write.csv2(res_all, paste("occurrences_cleaned/", name_data, "_all.csv", sep = ""), row.names = FALSE)  # Save to CSV
  
  # Generate species distribution plot (blue for AT_Database, red for GBIF)
  png(paste("graphics/", name_data, ".png", sep = ""))   # Save plot as PNG file
  plot(clim, ext = ext, legend = FALSE, main = name_data)  # Plot climate raster with species data
  points(res_w$X_round, res_w$Y_round, col = "blue", pch = 16)  # Plot AT_Database occurrences in blue
  points(res_g$X_round, res_g$Y_round, col = "red", pch = 16)  # Plot GBIF occurrences in red
  dev.off()  # Close the plot device
  
  # Create summary statistics for the data cleaning process
  single_stats <- data.frame(species = name_data, AT_DatabaseOrig = nrow(data_w), AT_DatabaseAggregated = nrow(res_w),
                             gbifOrig = nrow(data_g), gbifTime = nrow(data_time), gbifCertain = nrow(data_certain),
                             gbifCleaned = nrow(data_cl), gbifAggregated = nrow(res_g1), gbifAggregatedExtent = nrow(res_g),
                             AllCleaned = nrow(res_all))  
  all_stats <- rbind(all_stats, single_stats)  # Append statistics for the current species
}

# Save the final summary statistics to a CSV file
write.csv2(all_stats, "cleaning_stats1.csv", row.names = FALSE)  # Save the summary statistics to a CSV
