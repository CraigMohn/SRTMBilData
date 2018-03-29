library(tidyverse)
library(raster)
library(plotly)
library(htmlwidgets)
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(tigris)
library(sf)

source("C:/bda/SRTMBilData/R/functions.R")
source("C:/bda/SRTMBilData/R/regionDefs.R")
mapLibSelector <- 1  # 1=northAmerica+NE Pacific 1s, 2=Europe 3s, 3=Australia 3s
mapWindow <- NULL 
USStatevec <- "NV" # c("WA","OR") # c("MountainWest","CA","NM","AZ")
USParkvec <- NULL # <- c("CANY","CEBR","BRCA","ARCH") 
CAProvincevec <- NULL # "AB" #c("Maritimes","QC") #c("BC","AB","SK")
worldCountryvec <- NULL # <- c("DEU","AUT","CZE","CHE","FRA") #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html
showWater <- TRUE
waterLevel <- "area"  # "named","area","rivers","all"
showRoads <- TRUE
showCities <- TRUE

mapbuffer <- 1000 #5000     # meters, expands overall area
mapmergebuffer <- 200  # meters, expands categories before merging with others

cropbox <- raster::extent(-180, 170, -50, 60)  
#cropbox <- raster::extent(-180, 170, -50, 52.1) #  southern slice of BC,AB,SK 
res3dplot <- 3200
loadStateElevs <-  FALSE
writeElevFile <- TRUE
forceRes <- NULL
maxrastercells <- 500000000
highelevation <- 3000
vertscale <- 1.3
rglNAcolor <- "Blue"
drawRGL <- TRUE
saveRGL <- FALSE
drawPlotly <- FALSE
savePlotly <- FALSE
lon0to360=FALSE
workProj4 <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"

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
## mapWindow crops 
if (!is.null(mapWindow)) {
  mapcrop <- raster::extent(mapWindow)
  CP <- as(mapcrop, "SpatialPolygons")
  sp::proj4string(CP) <- workProj4
  mapcrop <- rgeos::gUnaryUnion(CP)
}
if (!is.null(USStatevec)) {
  USStatevec <- expandRegions(unique(toupper(USStatevec)),"US")  # US State abbrev all upper case
  statesInMap <- union(statesInMap,USStatevec)  
  mcrop <- rgeos::gUnaryUnion(sp::spTransform(tigris::counties(statesInMap),sp::CRS(workProj4))) 
  ## tigris returns a SpatialPolgonsDF and 
  ##    gUnaryUnion returns a SpatialPolygons
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
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
if (!is.null(CAProvincevec)) {
  CAProvincevec <- expandRegions(unique(toupper(CAProvincevec)),"CANADA")
  statesInMap <- union(statesInMap,CAProvincevec)
  canada <- raster::getData("GADM",country="CAN",level=1) # raster + spatial
  canada <-  sp::spTransform(canada,sp::CRS(workProj4))
  mcrop <- rgeos::gUnaryUnion(canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),])
  # simplify - BC coast is extremely complex
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1) 
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
mapcrop <- bufferUnion(mapcrop,mapbuffer=mapbuffer,NULL,simplifytol = 0)
CP <- as(cropbox, "SpatialPolygons")
sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))
mapcrop <- rgeos::gIntersection(mapcrop, CP, byid=TRUE)
plot(mapcrop)  #  which has CRS = workProj4

####################################################################################
if (loadStateElevs) {
  statesInMap <- USStatevec
  tmp <- loadStateElevData(USStatevec,CAProvincevec)  
  mapcrop <- raster::extent(cropbox)
  CP <- as(mapcrop, "SpatialPolygons")
  sp::proj4string(CP) <- workProj4
  mapcrop <- rgeos::gUnaryUnion(CP)
  spTown <- sp::spTransform(tmp[["spTown"]],sp::CRS(workProj4))
  spRoads <- sp::spTransform(tmp[["spRoads"]],sp::CRS(workProj4))
  spWaterA <- sp::spTransform(tmp[["spWaterA"]],sp::CRS(workProj4))
  spWaterL <- sp::spTransform(tmp[["spWaterL"]],sp::CRS(workProj4))
  m.sub <- tmp[["elevraster"]]
  if (!raster::compareCRS(raster::crs(m.sub),sp::CRS(workProj4)))
    m.sub <- raster::projectRaster(m.sub,crs=workProj4)
  
} else {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  if (!is.null(USStatevec)) {
    tmp <- USFeatures(USStatevec,workProj4,
                      showCities=showCities | writeElevFile,
                      showRoads=showRoads | writeElevFile,
                      showWater=showWater | writeElevFile)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  if (!is.null(CAProvincevec)) {
    tmp <- CAFeatures(CAProvincevec,workProj4,
                      showCities=showCities | writeElevFile,
                      showRoads=showRoads | writeElevFile,
                      showWater=showWater | writeElevFile)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  m.sub <- loadMapElevData(mapcrop)
  if (!raster::compareCRS(raster::crs(m.sub),sp::CRS(workProj4)))
    m.sub <- raster::projectRaster(m.sub,crs=workProj4)
}

if (lon0to360) m.sub <- raster::rotate(m.sub)
elevations <- m.sub

##  okay, elevs and shapefiles are set up

##################################################################################
if (writeElevFile & (length(statesInMap)==1)) {
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  if (nchunks == 1) {
    writeRaster(elevations,file=paste0(datadir,"/rasterfiles/",
                                       statesInMap[[1]],"elevs.grd"),
                overwrite=TRUE)   
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
    }
  }
  raster::shapefile(spRoads,filename=paste0(datadir,"/shapefiles/",
                                            statesInMap[[1]],"Roads.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterA,filename=paste0(datadir,"/shapefiles/",
                                             statesInMap[[1]],"WaterA.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterL,filename=paste0(datadir,"/shapefiles/",
                                             statesInMap[[1]],"WaterL.shp"),
                    overwrite=TRUE)
  raster::shapefile(spTown,filename=paste0(datadir,"/shapefiles/",
                                             statesInMap[[1]],"Town.shp"),
                    overwrite=TRUE)
}
#############################################################################


print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- max(1,floor(max(elevations@ncols,elevations@nrows)/res3dplot))
print(paste0("scaling raster down by a factor of ",sfact))
if (sfact > 1)
  elevations <- raster::aggregate(elevations,fact=sfact,fun=mean,
                           expand=TRUE,na.rm=FALSE)
print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))

if (drawRGL | drawPlotly) {

  yscale <- yRatio(elevations)
  mmmelev <- raster::as.matrix(elevations)
  if (drawPlotly) {
    mmm <- mmmelev[,ncol(mmmelev):1]  #  flip east/west since row 1 is top
    ax <- list(title="longitude",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    ay <- list(title="latitude",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    az <- list(title="elevation",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    p <- plotly::plot_ly(z = ~mmm,
                         colors = c("blue","yellow")) %>%
      plotly::add_surface(opacity=1.0) %>%
      plotly::layout(scene=list(xaxis=ax,yaxis=ay,zaxis=az,
                                aspectmode = "manual",
                                aspectratio = list(x=1,y=yscale,z=0.03*vertscale),
                                camera=list(up=c(0,1,0),
                                            eye=c(0,1.25,0)) ) )
    p
    if (savePlotly) htmlwidgets::saveWidget(p,paste0(mapoutputdir,"/",outputName," map.html"))
  }
  if (drawRGL) {  # mmmelev,cropbox, mapcrop.crs, showxx/spxx, palette, highelevation, yscale, vertscale
    x <- seq(1,length.out=nrow(mmmelev))
    y <- seq(1,length.out=ncol(mmmelev))
    
    CP <- as(doubleExtent(cropbox), "SpatialPolygons")
    sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))

    if (showCities | showWater | showRoads) {      
      print("building features")
      featurelayer <- raster::raster(nrow=nrow(elevations), ncol=ncol(elevations), 
                                     ext=extent(elevations), crs=crs(elevations))
      featurelayer[] <- 0
      if (showCities & !is.null(spTown)) {
        print("towns")
        spTown <- raster::crop(spTown, CP)
        featurelayer <- raster::rasterize(spTown,featurelayer,
                                          field=1,update=TRUE)
      }
      if (showWater) {
        spWaterL <- filterWaterL(spWaterL,level=waterLevel)
        if (!is.null(spWaterA)) {
          print("water polygons")
          spWaterA <- raster::crop(spWaterA, CP)
          featurelayer <- raster::rasterize(spWaterA,featurelayer,
                                            field=2,update=TRUE)
        }
        spWaterL <- filterWaterL(spWaterL,level=waterLevel)
        if (!is.null(spWaterL)) {
          print("water lines")
          spWaterL <- raster::crop(rgeos::gLineMerge(spWaterL), CP)
          featurelayer <- raster::rasterize(spWaterL,featurelayer,
                                            field=2,update=TRUE)
        }
      }
      if (showRoads & !is.null(spRoads)) {
        print("roads")
        spRoads <- raster::crop(rgeos::gLineMerge(spRoads), CP)
        featurelayer <- raster::rasterize(spRoads,featurelayer,
                                          field=3,update=TRUE)
      }
      features <- raster::as.matrix(featurelayer)
      print("features done")
    }
#   create a few palettes - PNW, islands, bright pastel
    terrcolors <- colorRampPalette(c("blue","turquoise","aquamarine",
                                     "palegreen","yellowgreen",
                                     "chartreuse","greenyellow","green",
                                     "limegreen","forestgreen","darkgreen",
                                     "yellow","gold","goldenrod",
                                     "sienna","brown","gray55",
                                     "gray65","gray80","gray90","white"))(206)
    terrcolors <- colorRampPalette(c("blue2",
                                     "palegreen","yellowgreen","lawngreen",
                                     "chartreuse","greenyellow","green",
                                     "limegreen","forestgreen","darkgreen",
                                     "yellow","gold","goldenrod",
                                     "sienna","brown","gray55",
                                     "gray65","gray80","gray90","white"))(206)
    plot(rep(1,206),col=terrcolors, pch=19,cex=2)
    
    tmpelev <- mmmelev/highelevation # don't worry about memory, not the constraint here
    tmpelev[is.na(tmpelev)] <- 0
    tmpelev <- sign(tmpelev)*sqrt(abs(tmpelev)) # f(0)=0, f(1)=1, f'(x>0) decreasing, reasonable for x<0
    tmpelev[mmmelev <= 0] <- 0  #  go off original 
    colidx <- floor(200*tmpelev) + 1
    colidx[colidx>201] <- 201
    colidx[colidx<1] <- 1
    colidx <- colidx + 5
    colidx[mmmelev <= 10]  <- 3
    col <- terrcolors[colidx]
    col[is.na(mmmelev)] <- gplots::col2hex(rglNAcolor)
    mmmelev[is.na(mmmelev)] <- -50
    if (showCities | showWater | showRoads) {
      col[features==1] <- gplots::col2hex("PeachPuff")
      col[features==2] <- gplots::col2hex("Blue")
      col[features==3] <- gplots::col2hex("Black")
    }
    rgl::par3d("windowRect"= c(100,100,1200,1000))
    userMatrix <- matrix(c(-0.02,-0.80,0.632,0,1,0,0.04,0,
                           -0.03,0.60,0.80,0,0,0,0,1),ncol=4,nrow=4)
    rgl::rgl.clear()
    rgl::surface3d(x,y,mmmelev,color=col)
    rgl::material3d(alpha=1.0,point_antialias=TRUE,smooth=TRUE,shininess=0)
    rgl::aspect3d(x=1,y=1/yscale,z=0.035*vertscale)
    rgl::rgl.clear("lights")
    rgl::rgl.light(theta = 0, phi = 15,
                   viewpoint.rel=TRUE, specular="black")
    rgl::rgl.viewpoint(userMatrix=userMatrix,type="modelviewpoint")
    pan3d(2)  # right button for panning, doesn't play well with zoom)
    if (saveRGL) rgl::writeWebGL(dir=paste0(mapoutputdir), filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))
  }
}
