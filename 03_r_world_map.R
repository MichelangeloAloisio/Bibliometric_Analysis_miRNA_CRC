library(classInt)
library(RColorBrewer)
library(rworldmap)

# Load and merge data
data("countryExData", envir = environment(), package = "rworldmap")
country_data <- read.csv(
  "D:/Desktop/progetti_v1/20_Serino_colon_bibliometrix/01_Analisi_Bibliometrica/01_extract_miRNA/Country_Production_bibliometrix.csv", 
  sep = ";", 
  header = TRUE
)

# Merge data with map
spdf <- joinCountryData2Map(
  country_data, 
  joinCode = "NAME", 
  nameJoinColumn = "Country", 
  mapResolution = 'coarse', 
  verbose = TRUE
)

# Verify column names (important!)
print(names(spdf))  # Confirm "Frequency.Percentage" exists

# Filter out NA values (if any)
valid_data <- spdf$Frequency.Percentage[!is.na(spdf$Frequency.Percentage)]

# Create quantile-based classification intervals
classInt <- classIntervals(valid_data, n = 4, style = "quantile")

# Extract breaks for legend
catMethod <- classInt$brks

# Define color palette
basePalette <- colorRampPalette(c("#FFEDCC", "#F8C471", "#F1948A", "#E74C3C"))(4)
colourPalette <- adjustcolor(basePalette, alpha.f = 0.85)

# Plot the map using Frequency.Percentage
mapParams <- mapCountryData(
  spdf, 
  nameColumnToPlot = "Frequency.Percentage", 
  addLegend = FALSE, 
  catMethod = catMethod, 
  colourPalette = colourPalette
)

legend("bottomright", 
       legend = paste0(
         "[", 
         round(catMethod[-length(catMethod)], 3),  # Rounded to 3 decimals
         " - ", 
         round(catMethod[-1], 3),                  # Rounded to 3 decimals
         "]%"
       ),
       fill = colourPalette,
       border = "black",
       title = "Publication Percentage",
       cex = 0.8,
       bg = "white")

