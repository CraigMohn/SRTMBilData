library(tidyverse)
library(raster)
library(htmlwidgets)
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(tigris)
library(sf)
library(rmapshaper)

source("C:/bda/SRTMBilData/R/functions.R")
source("C:/bda/SRTMBilData/R/featurefunctions.R")
source("C:/bda/SRTMBilData/R/displayfunctions.R")
source("C:/bda/SRTMBilData/R/regionDefs.R")
mapWindow <- USStatevec <- USParkvec <- CAProvincevec <- worldCountryvec <- NULL 
mapLibSelector <- 1  # 1=northAmerica+NE Pacific 1s, 2=Europe 3s, 3=Australia 3s

#mapWindow <- c(-122.5,-121.9,37.6,38.1)     # East Bay 
USStatevec <- c("WA") #  c("MountainWest","CA","NM","AZ")
#USParkvec <- c("CANY","CEBR","BRCA","ARCH") 
#CAProvincevec <- "BC" # c("BC","AB","SK")
#worldCountryvec <-  c("DEU","AUT","CZE","CHE","FRA") #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html
showWater <- TRUE
waterLevel <- "riversplus" # "named","lakes",area","rivers","riversplus","all"
showRoads <- TRUE
showCities <- TRUE

mapbuffer <- 5000     # meters, expands overall area
mapmergebuffer <- 0 #200  # meters, expands categories before merging with others

cropbox <- raster::extent(-180, 170, -50, 60)
if (!is.null(mapWindow)) cropbox <- raster::extent(mapWindow)
#cropbox <- raster::extent(-160.25, -154.8, 18.9, 22.25) # hawaii main islands only
#cropbox <- raster::extent(-180, 170, -50, 52.1) #  southern slice of BC,AB,SK 
res3dplot <- 3000
loadStateElevs <-  FALSE
writeElevFile <- TRUE
loadStateFeatureRaster <- FALSE
writeStateFeatureRaster <- TRUE
forceRes <- NULL
maxrastercells <- 500000000
highelevation <- 3000
vertscale <- 2.25   
rglColorScheme <- "default"
rglNAcolor <- "Blue"
rglNegcolor <- "Red"
citycolor <- "white"
drawRGL <- TRUE
saveRGL <- TRUE
lon0to360=FALSE
workProj4 <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"
options(tigris_use_cache = TRUE)

datadir <- "c:/bda"                    #  base dir - subs include shapefiles, rasterfile 
mapoutputdir <- "c:/bda/maps3d"        #  map file output location
NAmericaDataDir <- "c:/bda/NorthAmerica"   #  zip file input subdirectories location
EuropeDataDir <- "c:/bda/Europe 3s"        #  zip file input subdirectories location
AustraliaDataDir <- "c:/bda/Australia 3s"        #  zip file input subdirectories location
mapDataDir <- c(NAmericaDataDir,EuropeDataDir,AustraliaDataDir)[[mapLibSelector]]
resstr <- c("_1arc_v3_bil","_3arc_v2_bil","_1arc_v3_bil")[[mapLibSelector]]

#################################################################################
assign("last.warning", NULL, envir = baseenv())
if (loadStateElevs & writeElevFile) stop("don't ask me to read the data and then write it")
outputName <- paste0(c(USStatevec,CAProvincevec,USParkvec),collapse="-")

####    set up cropping mask file and vector of states/etc
mapcrop <- NULL
statesInMap <- NULL
if (!is.null(USStatevec)) {
  USStatevec <- expandRegions(unique(toupper(USStatevec)),"US")  # US State abbrev all upper case
  statesInMap <- union(statesInMap,USStatevec)  
  mcrop <- tigris::counties(statesInMap) %>% 
    rgeos::gUnaryUnion(.) %>% 
    sp::spTransform(.,sp::CRS(workProj4))
  ## tigris returns a SpatialPolgonsDF and gUnaryUnion returns a SpatialPolygons
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
}
if (!is.null(CAProvincevec)) {
  CAProvincevec <- expandRegions(unique(toupper(CAProvincevec)),"CANADA")
  statesInMap <- union(statesInMap,CAProvincevec)
  canada <- raster::getData("GADM",country="CAN",level=1) %>%
    sp::spTransform(.,sp::CRS(workProj4))
  mcrop <- canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),] %>%
    rgeos::gUnaryUnion(.)
  # simplify - BC coast is extremely complex
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1) 
}
if (!is.null(USParkvec)) {
  USParkvec <- unique(toupper(USParkvec))
  statesInMap <- union(statesInMap,USParkvec)
  pdir <- "nps_boundary"
  pfile <- paste0(pdir,".shp")
  parkareas <- sf::st_read(paste0(datadir,"/",pdir,"/",pfile))  # sf dataframe
  #parkareas <- rgeos::gUnaryUnion(parkareas)
  #  parknames <- parkareas[,c("UNIT_NAME","UNIT_CODE")]
  #  sf::st_geometry(parknames) <- NULL
  #  parknames <- parknames[order(parknames$UNIT_CODE),] 
  #  write.csv(parknames,paste0(datadir,"/parknames.csv"))
  parkareas <-  sp::spTransform(parkareas[parkareas$UNIT_CODE %in% USParkvec,"geometry"],
                                sp::CRS(workProj4)) 
  mcrop <- as(parkareas, "Spatial")  
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
}
if (!is.null(worldCountryvec)) { 
  worldCountryvec <- unique(toupper(worldCountryvec))
  statesInMap <- union(statesInMap,worldCountryvec)
  mcrop <- NULL
  for (c in worldCountryvec) {
    cmap <- raster::getData("GADM",country=worldCountryvec,level=0) # raster + spatial
    cmap <- rgeos::gUnaryUnion(cmap) 
    cmap <-  sp::spTransform(cmap,sp::CRS(workProj4))
    if (is.null(mcrop)) {
      mcrop <- cmap
    } else {
      mcrop<- rgeos::gUnaryUnion(raster::union(mcrop,cmap))
    }
  }
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1)
}
## mapWindow - overwrite the map used for cropping, 
##     but not the list of states/provinces to load
if (!is.null(mapWindow)) {
  mapcrop <- raster::extent(mapWindow)
  CP <- as(mapcrop, "SpatialPolygons")
  sp::proj4string(CP) <- workProj4
  mapcrop <- rgeos::gUnaryUnion(CP)
}
mapcrop <- bufferUnion(mapcrop,mapbuffer=mapbuffer,NULL,simplifytol = 0)

#   now crop by the cropbox 
if (writeElevFile) warning("cropping map when saving raster/shapefiles.")
CP <- as(cropbox, "SpatialPolygons")
sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))
mapcrop <- rgeos::gIntersection(mapcrop, CP, byid=TRUE)
plot(mapcrop)  #  which has CRS = workProj4

####################################################################################
if (loadStateElevs) {
  statesInMap <- USStatevec
  print(system.time(tmp <- loadStateElevData(USStatevec,CAProvincevec))[[3]])
  spTown <- sp::spTransform(tmp[["spTown"]],sp::CRS(workProj4))
  spRoads <- sp::spTransform(tmp[["spRoads"]],sp::CRS(workProj4))
  spWaterA <- sp::spTransform(tmp[["spWaterA"]],sp::CRS(workProj4))
  spWaterL <- sp::spTransform(tmp[["spWaterL"]],sp::CRS(workProj4))
  elevations <- tmp[["elevraster"]]
  if (!raster::compareCRS(raster::crs(elevations),sp::CRS(workProj4)))
    elevations <- raster::projectRaster(elevations,crs=workProj4)
  print(system.time(elevations <- raster::mask(elevations,mapcrop))[[3]])
} else {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  if (!is.null(USStatevec)) {
    tmp <- USFeatures(USStatevec,workProj4,
                      writeShapefiles=writeElevFile,
                      shapefiledir=paste0(datadir,"/shapefiles/"))
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  if (!is.null(CAProvincevec)) {
    tmp <- CAFeatures(CAProvincevec,workProj4)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  elevations <- loadMapElevData(mapcrop)
  if (!raster::compareCRS(raster::crs(elevations),sp::CRS(workProj4)))
    elevations <- raster::projectRaster(elevations,crs=workProj4)
}
spTown <- raster::intersect(spTown,mapcrop)

if (writeStateFeatureRaster) 
  featureStack <- buildFeatureStack(elevations,mapshape=mapcrop,
                                    spTown,spWaterA,spWaterL,spRoads,
                                    fastAreas=FALSE,
                                    polySimplify=0,polyMethod="vis",
                                    polyWeighting=0.85,polySnapInt=0.0001) %>%
    raster::mask(.,mapcrop)
  
##  okay, elevs and shapefiles are set up

##################################################################################
if (writeElevFile & (length(statesInMap)==1)) {
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  print(paste0("saving raster data in ",nchunks," slices"))
  if (nchunks == 1) {
    writeRaster(elevations,file=paste0(datadir,"/rasterfiles/",
                                       statesInMap[[1]],"elevs.grd"),
                overwrite=TRUE)   
    if (writeStateFeatureRaster) 
      writeRaster(featureStack,file=paste0(datadir,"/rasterfiles/",
                  statesInMap[[1]],"features.grd"),
                  datatype='INT1S',overwrite=TRUE)
  } else {
    nrowchunk <- ceiling(raster::nrow(elevations)/nchunks)
    for (chunk in 1:nchunks) {
      chunkcrop <- raster::extent(raster::xmin(elevations),
                                  raster::xmax(elevations),
                                  raster::ymin(elevations) + 
                                    raster::yres(elevations)*
                                    ((chunk-1)*nrowchunk - 10),
                                  raster::ymin(elevations) + 
                                    raster::yres(elevations)*
                                    (chunk*nrowchunk - 10)  )
      chunkcrop <- as(chunkcrop, "SpatialPolygons")
      sp::proj4string(chunkcrop) <- sp::proj4string(elevations)
      chunkraster <- raster::trim(raster::crop(elevations,chunkcrop))
      writeRaster(chunkraster,
                  file=paste0(datadir,"/rasterfiles/",
                              statesInMap[[1]],"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  overwrite=TRUE)    
      if (writeStateFeatureRaster) {
        chunkfeature <- raster::trim(raster::crop(featureStack,chunkcrop))
        writeRaster(chunkfeature,
                    file=paste0(datadir,"/rasterfiles/",
                                statesInMap[[1]],"features",
                                stringr::str_pad(chunk,2,pad="0"),".grd"),
                    datatype='INT1S',overwrite=TRUE)    
      }
    }
  }
}
#############################################################################
# crop raster after write
elevations <- raster::crop(elevations,mapcrop)
if (lon0to360) elevations <- raster::rotate(elevations)

print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- ceiling(sqrt(elevations@ncols*elevations@nrows/(res3dplot*res3dplot)))
print(paste0("scaling raster down by a factor of ",sfact))
if (sfact > 1)
  print(system.time(
    elevations <- raster::aggregate(elevations,fact=sfact,fun=mean,
                                    expand=TRUE,na.rm=FALSE)
  )[[3]])
  print(system.time(
    featureStack <- raster::aggregate(featureStack,fact=sfact,fun=max,
                                  expand=TRUE,na.rm=TRUE)
  )[[3]])
print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))

if (drawRGL) {
  pspTown <- NULL
  pspRoads <- NULL
  pspWaterA <- NULL
  pspWaterL <- NULL
  if (showCities) pspTown <- spTown
  if (showRoads) pspRoads <- spRoads
  if (showWater) {
    pspWaterA <- spWaterA
    pspWaterL <- spWaterL
  }
  print(system.time(
  draw3DMapRgl(elevations,mapPolygon=mapcrop,
               spTown=pspTown,spRoads=pspRoads,
               spWaterA=pspWaterA,spWaterL=pspWaterL,
               waterLevel,colors=rglColorScheme,
               maxElev=highelevation,vScale=vertscale/(1.0 + log(sfact)),
               rglNAcolor=rglNAcolor,
               rglNegcolor=rglNegcolor,
               citycolor=citycolor,
               saveRGL=saveRGL,
               mapoutputdir=mapoutputdir,outputName=outputName,
               fastAreas=FALSE,polySimplify=0.0,
               polyMethod="vis",polyWeighting=0.88,
               polySnapInt=0.0001) 
)[[3]])
}
