library(classInt)
library(RColorBrewer)
library(rworldmap)

# Carica i dati e uniscili alla mappa
data("countryExData", envir=environment(), package="rworldmap")
country_data <- read.csv("D:/Desktop/progetti_v1/20_Serino_colon_bibliometrix/01_Analisi_Bibliometrica/01_extract_miRNA/Country_Production_bibliometrix.csv", sep=';', header=TRUE)
spdf <- joinCountryData2Map(country_data, joinCode="NAME", nameJoinColumn="Country", mapResolution='coarse', verbose=TRUE)

# Controlla eventuali paesi non corrispondenti e statistiche Freq
summary(spdf)
summary(spdf$Freq)

# Crea intervalli di classificazione (quantili)
classInt <- classIntervals(spdf[["Freq"]], n=4, style="quantile")
catMethod <- classInt[["brks"]]

# Definisci la palette di colori pastello: dal giallo chiaro al rosso chiaro
basePalette <- colorRampPalette(c("#FFEDCC", "#F8C471", "#F1948A", "#E74C3C"))(4) # Palette pastello

# Rendi i colori leggermente più opachi
colourPalette <- adjustcolor(basePalette, alpha.f = 0.85) # 


# Mappa i dati
mapParams <- mapCountryData(spdf, nameColumnToPlot="Freq", addLegend=FALSE, catMethod=catMethod, colourPalette=colourPalette)

# Aggiungi la legenda manualmente
legend("bottomright", # Posizione della legenda (puoi modificare)
       legend = paste0("[", round(catMethod[-length(catMethod)], 2), " - ", round(catMethod[-1], 2), "]"),
       fill = colourPalette,
       border = "black",
       title = "Frequency Intervals",
       cex = 0.8) # Dimensione del testo





