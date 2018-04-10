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
mapWindow <- USStatevec <- USParkvec <- CAProvincevec <- worldCountryvec <- NULL 
mapLibSelector <- 1  # 1=northAmerica+NE Pacific 1s, 2=Europe 3s, 3=Australia 3s

mapWindow <- c(-122.7,-121.8,37.4,38.3)     # East Bay
USStatevec <- c("CA") #  c("MountainWest","CA","NM","AZ")
#USParkvec <- c("CANY","CEBR","BRCA","ARCH") 
#CAProvincevec <- "BC" # c("BC","AB","SK")
#worldCountryvec <-  c("DEU","AUT","CZE","CHE","FRA") #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html
showWater <- TRUE
waterLevel <- "riversplus"  # "named","lakes",area","rivers","riversplus","all"
showRoads <- TRUE
showCities <- TRUE

mapbuffer <- 1000 #5000     # meters, expands overall area
mapmergebuffer <- 200  # meters, expands categories before merging with others

cropbox <- raster::extent(-180, 170, -50, 60)  
if (!is.null(mapWindow)) cropbox <- raster::extent(mapWindow)
#cropbox <- raster::extent(-160.25, -154.8, 18.9, 22.25) # hawaii main islands only
#cropbox <- raster::extent(-180, 170, -50, 52.1) #  southern slice of BC,AB,SK 
res3dplot <- 3200
loadStateElevs <-  FALSE
writeElevFile <- FALSE
forceRes <- NULL
maxrastercells <- 500000000
highelevation <- 3000
vertscale <- 1.8   
rglNAcolor <- "Blue"
citycolor <- "White"
drawRGL <- TRUE
saveRGL <- TRUE
drawPlotly <- FALSE
savePlotly <- FALSE
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

#   now crop by the cropbox unless saving
if (!writeElevFile) {
  CP <- as(cropbox, "SpatialPolygons")
  sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))
  mapcrop <- rgeos::gIntersection(mapcrop, CP, byid=TRUE)
}
plot(mapcrop)  #  which has CRS = workProj4

####################################################################################
if (loadStateElevs) {
  statesInMap <- USStatevec
  print(system.time(tmp <- loadStateElevData(USStatevec,CAProvincevec))[[3]])
  spTown <- sp::spTransform(tmp[["spTown"]],sp::CRS(workProj4))
  spRoads <- sp::spTransform(tmp[["spRoads"]],sp::CRS(workProj4))
  spWaterA <- sp::spTransform(tmp[["spWaterA"]],sp::CRS(workProj4))
  spWaterL <- sp::spTransform(tmp[["spWaterL"]],sp::CRS(workProj4))
  m.sub <- tmp[["elevraster"]]
  if (!raster::compareCRS(raster::crs(m.sub),sp::CRS(workProj4)))
    m.sub <- raster::projectRaster(m.sub,crs=workProj4)
  m.sub <- raster::mask(m.sub,mapcrop)
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
  m.sub <- loadMapElevData(mapcrop)
  if (!raster::compareCRS(raster::crs(m.sub),sp::CRS(workProj4)))
    m.sub <- raster::projectRaster(m.sub,crs=workProj4)
}
spTown <- sxdfMask(spTown,mapcrop)
elevations <- m.sub


##  okay, elevs and shapefiles are set up

##################################################################################
if (writeElevFile & (length(statesInMap)==1)) {
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  print(paste0("saving raster data in ",nchunks," slices"))
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
}
#############################################################################
# crop raster after write
elevations <- raster::crop(elevations,mapcrop)
if (lon0to360) elevations <- raster::rotate(elevations)

print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- max(1,floor(max(elevations@ncols,elevations@nrows)/res3dplot))
print(paste0("scaling raster down by a factor of ",sfact))
if (sfact > 1)
  print(system.time(
    elevations <- raster::aggregate(elevations,fact=sfact,fun=mean,
                                    expand=TRUE,na.rm=FALSE)
  )[[3]])
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
    if (savePlotly) 
      htmlwidgets::saveWidget(p,paste0(mapoutputdir,"/",outputName," map.html"))
  }
  if (drawRGL) {  # mmmelev,mapcrop, showxx/spxx, palette, highelevation, yscale, vertscale
    x <- seq(1,length.out=nrow(mmmelev))
    y <- seq(1,length.out=ncol(mmmelev))
    
    if (showCities | showWater | showRoads) {      
      CP <- as(doubleExtent(elevations), "SpatialPolygons")
      sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))

      print("building features")
      featurelayer <- raster::raster(nrow=nrow(elevations),ncol=ncol(elevations), 
                                     ext=extent(elevations),crs=crs(elevations))
      featurelayer[] <- 0
      print(featurelayer@data@min)
      print(featurelayer@data@max)
      if (showCities & !is.null(spTown)) {
        print("towns")
        tspTown <- raster::intersect(spTown, mapcrop)
        if (!is.null(tspTown)) {
          plot(tspTown,col=citycolor)
          for (i in 1:nrow(tspTown)) {
            print(tspTown@data[i,"NAME"])
            plot(tspTown[i,],col=citycolor)
            featurelayer <- raster::rasterize(tspTown[i,],featurelayer,
                                              field=1,update=TRUE)
          }
        }
        print(featurelayer@data@min)
        print(featurelayer@data@max)
      }
      if (showWater & !is.null(spWaterA)) {
        tspWaterA <- filterWaterA(spWaterA,level=waterLevel) # need new copy since not in function
        tspWaterA <- raster::intersect(tspWaterA, mapcrop)
        if (!is.null(tspWaterA)) {
          print("water polygons")
          plot(tspWaterA, col="blue")
          featurelayer <- raster::rasterize(tspWaterA,featurelayer,
                                          field=2,update=TRUE)
        }
        print(featurelayer@data@min)
        print(featurelayer@data@max)
      }
      if (showWater & !is.null(spWaterL)) {
        tspWaterL <- filterWaterL(spWaterL,level=waterLevel)  # need new copy since not in function
        if (!is.null(tspWaterL)) {
          print("water lines")
          tspWaterL <- raster::crop(rgeos::gLineMerge(tspWaterL), CP)
          plot(tspWaterL, col="blue")
          featurelayer <- raster::rasterize(tspWaterL,featurelayer,
                                            field=2,update=TRUE)
        }
        print(featurelayer@data@min)
        print(featurelayer@data@max)
      }
      if (showRoads & !is.null(spRoads)) {
        print("roads")
        tspRoads <- raster::crop(rgeos::gLineMerge(spRoads), CP)
        plot(tspRoads)
        featurelayer <- raster::rasterize(tspRoads,featurelayer,
                                          field=3,update=TRUE)
        print(featurelayer@data@min)
        print(featurelayer@data@max)
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
    ### standard
    terrcolors <- colorRampPalette(c("blue","darkturquoise","turquoise","aquamarine",
                                     "palegreen","greenyellow","lawngreen",
                                     "chartreuse","green","springgreen",
                                     "limegreen","forestgreen","darkgreen",
                                     "olivedrab","darkkhaki","darkgoldenrod",
                                     "sienna","brown","saddlebrown","rosybrown",
                                     "gray35","gray45","gray55",
                                     "gray65","gray70","gray75","gray85"))(206)
    ### beaches
#    terrcolors <- colorRampPalette(c("blue","bisque1","bisque2","bisque3",
#                                     "palegreen","greenyellow","lawngreen",
#                                     "chartreuse","green","springgreen",
#                                     "limegreen","forestgreen","darkgreen",
#                                     "olivedrab","darkkhaki","darkgoldenrod",
#                                     "sienna","brown","saddlebrown","rosybrown",
#                                     "gray35","gray45","gray55",
#                                     "gray65","gray70","gray75","gray85"))(206)
    plot(rep(1,206),col=terrcolors, pch=19,cex=2)
    
    tmpelev <- mmmelev/highelevation # don't worry about memory, not the constraint here
    tmpelev[is.na(tmpelev)] <- 0
    tmpelev <- sign(tmpelev)*sqrt(abs(tmpelev)) # f(0)=0, f(1)=1, f'(x>0) decreasing, reasonable for x<0
    tmpelev[mmmelev <= 0] <- 0  #  go off original 
    colidx <- floor(200*tmpelev) + 1
    colidx[colidx>201] <- 201
    colidx[colidx<1] <- 1
    colidx <- colidx + 5
    colidx[mmmelev == 0]  <- 1
    col <- terrcolors[colidx]
    col[is.na(mmmelev)] <- gplots::col2hex(rglNAcolor)
    mmmelev[is.na(mmmelev)] <- -50
    if (showCities | showWater | showRoads) {
      col[features==1] <- gplots::col2hex(citycolor)
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
    if (saveRGL) 
      rgl::writeWebGL(dir=paste0(mapoutputdir), 
                      filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))
  }
}
