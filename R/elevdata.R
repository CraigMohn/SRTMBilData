library(tidyverse)
library(raster)
library(htmlwidgets)
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(tigris)
library(sf)
library(velox)
library(rmapshaper)
library(viridis)

source("C:/bda/SRTMBilData/R/drawMapRGL.R")
source("C:/bda/SRTMBilData/R/functions.R")
source("C:/bda/SRTMBilData/R/featurefunctions.R")
source("C:/bda/SRTMBilData/R/displayfunctions.R")
source("C:/bda/SRTMBilData/R/mapMask.R")
source("C:/bda/SRTMBilData/R/readwrite.R")
source("C:/bda/SRTMBilData/R/readSRTMdata.R")
source("C:/bda/SRTMBilData/R/regionDefs.R")

mapLibSelector <- 1  # 1=northAmerica+NE Pacific 1s, 2=Europe 3s, 3=Australia 3s

#mapWindow <- c(-122.5,-121.9,37.6,38.1)     # East Bay 
#mapWindow <- c(-122.7,-121.8,37.0,38.3)     # SF Bay
#USStatevec <- c("WA") #  c("MountainWest","CA","NM","AZ")
#USParkvec <- c("CANY","CEBR","BRCA","ARCH") 
#CAProvincevec <- "BC" # c("BC","AB","SK")
#worldCountryvec <-  c("DEU","AUT","CZE","CHE","FRA") #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html

#cropbox <- raster::extent(-180, 170, -50, 60)
#if (!is.null(mapWindow)) cropbox <- raster::extent(mapWindow)
#cropbox <- raster::extent(-160.25, -154.8, 18.9, 22.25) # hawaii main islands only
#cropbox <- raster::extent(-180, 170, -50, 52.1) #  southern slice of BC,AB,SK 
options(tigris_use_cache = TRUE)

datadir <- "c:/bda"                        #  base dir - subs include shapefiles, rasterfile 
parkDir <- paste0(datadir,"/nps_boundary")
rasterDir <- paste0(datadir,"/rasterfiles")
shapefileDir <- paste0(datadir,"/shapefiles")
mapoutputdir <- "c:/bda/maps3d"            #  map file output location
NAmericaDataDir <- "c:/bda/NorthAmerica"   #  zip file input subdirectories location
EuropeDataDir <- "c:/bda/Europe 3s"        #  zip file input subdirectories location
AustraliaDataDir <- "c:/bda/Australia 3s"        #  zip file input subdirectories location
mapDataDir <- c(NAmericaDataDir,EuropeDataDir,AustraliaDataDir)[[mapLibSelector]]
resstr <- c("_1arc_v3_bil","_3arc_v2_bil","_1arc_v3_bil")[[mapLibSelector]]

#################################################################################



drawMapRGL(USStatevec="OR",
           elevDataSource="SRTM",
           featureDataSource="Shapefiles",
           writeElevFile=TRUE,writeFeatureFile=TRUE,
           rasterDir=rasterDir,mapDataDir=mapDataDir,
           shapefileDir=shapefileDir,
           saveRGL=TRUE,mapoutputdir=mapoutputdir,
           outputName="OR")



