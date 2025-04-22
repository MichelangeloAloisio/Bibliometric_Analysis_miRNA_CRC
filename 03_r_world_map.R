library(classInt)
library(RColorBrewer)
library(rworldmap)

# Load and join data (with minor fixes)
data("countryExData", envir = environment(), package = "rworldmap")
country_data <- read.csv("D:/Desktop/progetti_v1/20_Serino_colon_bibliometrix/01_Analisi_Bibliometrica/01_extract_miRNA/Country_Production_bibliometrix.csv", 
                        sep = ';', header = TRUE)

spdf <- joinCountryData2Map(country_data, 
                           joinCode = "NAME", 
                           nameJoinColumn = "Country", 
                           mapResolution = 'coarse', 
                           verbose = TRUE)

# Handle missing values (if any)
spdf <- spdf[!is.na(spdf$Freq), ]

# Create class intervals (4 quantiles)
classInt <- classIntervals(spdf[["Freq"]], n = 4, style = "quantile")
catMethod <- classInt$brks  # Extract breaks as a numeric vector

# Choose a color palette (4 colors for 4 classes)
colourPalette <- brewer.pal(4, "YlOrRd")  # Use "YlOrRd" for a clear gradient

# Plot the map
mapParams <- mapCountryData(
  spdf,
  nameColumnToPlot = "Freq",
  addLegend = TRUE,
  catMethod = catMethod,
  colourPalette = colourPalette
)

# Optional: Customize legend (if needed)
do.call(addMapLegend, c(mapParams, 
                       legendLabels = "all", 
                       legendWidth = 0.5, 
                       legendIntervals = "data", 
                       legendMar = 8))
