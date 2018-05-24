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
library(openStreetMap)
library(viridis)

source("C:/bda/SRTMBilData/R/drawMapRGL.R")
source("C:/bda/SRTMBilData/R/featuresForElevations.R")
source("C:/bda/SRTMBilData/R/functions.R")
source("C:/bda/SRTMBilData/R/featurefunctions.R")
source("C:/bda/SRTMBilData/R/displayfunctions.R")
source("C:/bda/SRTMBilData/R/mapMask.R")
source("C:/bda/SRTMBilData/R/readwrite.R")
source("C:/bda/SRTMBilData/R/readSRTMdata.R")
source("C:/bda/SRTMBilData/R/regionDefs.R")

mapLibSelector <- 1  # 1=northAmerica+NE Pacific 1s, 2=Europe 3s, 3=Australia 1s

#USStatevec <- c("WA") #  c("MountainWest","CA","NM","AZ")
#USParkvec <- c("CANY","CEBR","BRCA","ARCH") 
#CAProvincevec <- "BC" # c("BC","AB","SK")
#worldCountryvec <-  c("DEU","AUT","CZE","CHE","FRA") #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html

#cropbox <- raster::extent(-180, 170, -50, 60)
#if (!is.null(mapWindow)) cropbox <- raster::extent(mapWindow)
#cropbox <- raster::extent(-160.25, -154.8, 18.9, 22.25) # hawaii main islands only
cropbox <- raster::extent(-180, 170, -50, 54.7) #  southern slice of BC,AB,SK 
options(tigris_use_cache = TRUE)

datadir <- "c:/bda"                        #  base dir - subs include shapefiles, rasterfile 
parkDir <- paste0(datadir,"/nps_boundary")
rasterDir <- paste0(datadir,"/rasterfiles")
shapefileDir <- paste0(datadir,"/shapefiles")
mapoutputdir <- "c:/bda/maps3d"            #  map file output location
NAmericaDataDir <- "c:/bda/NorthAmerica"   #  zip file input subdirectories location
EuropeDataDir <- "c:/bda/Europe 3s"        #  zip file input subdirectories location
AustraliaDataDir <- "c:/bda/Australia 1s"  #  zip file input subdirectories location
mapDataDir <- c(NAmericaDataDir,EuropeDataDir,AustraliaDataDir)[[mapLibSelector]]
resstr <- c("_1arc_v3_bil","_3arc_v2_bil","_1arc_v1_bil")[[mapLibSelector]]

#################################################################################

mapWindow <- NULL
# mapWindow <- c(-156.80,-155.97,20.50,21.05) # Maui
# mapWindow <- c(-159.90,-159.15,21.75,22.35) # Kauai 
# mapWindow <- c(-156.10,-154.75,18.85,20.30) # Big Island
# mapWindow <- c(-81.4,-80.0,36.8,37.6)       # Giles Cty/Blacksburg Area 
# mapWindow <- c(-123.25,-121.5,46.75,48.1)   # Seattle Area 
# mapWindow <- c(-123.2,-122.4,48.3,48.8)     # San Juans
# mapWindow <- c(-122.2,-121.7,47.4,47.8)     # Samm Area 

drawMapRGL(USStatevec=NULL, 
           CAProvincevec="NS",
           mapWindow=mapWindow,
           #cropbox=cropbox,
           elevDataSource="SRTM",
           featureDataSource="Shapefiles",
           townLevel=5,roadLevel=3,waterALevel=4,waterLLevel=5,
           vScale=1.25,maxElev=2000,
           rglColorScheme="default",
           #citycolor="purple",
           writeElevFile=TRUE,
           writeFeatureFile=TRUE,
           rasterDir=rasterDir,mapDataDir=mapDataDir,
           shapefileDir=shapefileDir,includeAllRoads=TRUE,
           maxrastercells=200000000,saveRGL=FALSE,res3d=2800,
           maxRasterize=100000,
           mapoutputdir=mapoutputdir,
           #rasterFileSetWriteNames="SanJuans",
           outputName="NS")

stop()
elevations <- elevationsToRaster(rasterFileSetName="SK",
                               USStatevec=NULL,CAProvincevec="SK",
                               cropbox=cropbox,
                               rasterDir=rasterDir,mapDataDir=mapDataDir,
                               mapbuffer=0,mapmergebuffer=0,
                               maxrastercells=250000000,
                               workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                               resstr="_1arc_v3_bil") 
    

featuresForElevations(rasterFileSetName="NS",
                      rasterDir=rasterDir,
                      shapefileDir=shapefileDir,
                      USStatevec="NS",
                      CAProvincevec=NULL,
                      featureDataSource="Shapefiles",
                      includeAllRoads=TRUE,
                      workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                      maxRasterize=100000,
                      polySimplify=0.0,polyMethod="vis", 
                      polyWeighting=0.85,polySnapInt=0.0001) 

drawMapRGL(USStatevec=NULL, 
           CAProvincevec=c("NS"),
           mapWindow=mapWindow,
           #cropbox=cropbox,
           elevDataSource="Raster",
           featureDataSource="Raster",
           townLevel=5,roadLevel=4,waterALevel=4,waterLLevel=5,
           vScale=1.25,maxElev=2000,
           rglColorScheme="default",
           rasterDir=rasterDir,
           saveRGL=TRUE,res3d=2800,
           mapoutputdir=mapoutputdir,
           outputName="NS")

  