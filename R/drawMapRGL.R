drawMapRGL <- function(USStatevec=NULL,CAProvincevec=NULL,USParkvec=NULL,
                       worldCountryvec=NULL,mapWindow=NULL,
                       routeSL=NULL,cropbox=NULL,
                       res3dplot=2500,maxElev=3000,vScale=1.5,
                       townLevel=3,roadLevel=4,waterALevel=4,waterLLevel=5,
                       rglColorScheme="default",
                       rglNAcolor="Blue",rglNegcolor="Red",
                       citycolor="SlateGray",watercolor="Blue",
                       roadcolor="Black",glaciercolor="White",
                       drawRGL=TRUE,
                       saveRGL=FALSE,mapoutputdir=NULL,outputName=NULL,
                       workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                       elevDataSource="Raster",useRasterFileSets=NULL,
                       rasterFileSetNames=NULL,rasterFileSetWriteNames=NULL,
                       featureDataSource="Shapefiles",
                       writeElevFile=FALSE,writeFeatureFile=FALSE,
                       writeShapefiles=TRUE,includeAllRoads=FALSE,
                       rasterDir=NULL,mapDataDir=NULL,shapefileDir=NULL,parkdir=NULL,
                       resstr="_1arc_v3_bil",
                       mapbuffer=0,mapmergebuffer=0,
                       maxrastercells=200000000,
                       polySimplify=0.0,polyMethod="vis", 
                       polyWeighting=0.85,polySnapInt=0.0001) {

  if (elevDataSource=="Raster" & writeElevFile &
           is.null(rasterFileSetWriteNames)) {
    warning("will not overwrite elev rasterfile that was source")
    writeElevFile <- FALSE
  }
  if (featureDataSource=="Raster" & writeFeatureFile &
           is.null(rasterFileSetWriteNames)){
    warning("will not overwrite feature rasterfile that was source")
    writeFeatureFile <- FALSE
  }

  featureFilter <- c(townLevel,roadLevel,waterALevel,waterLLevel)
  
  ##  to do:  add list of saved elev/feature rasters
  ##          default mapWindow to user-supplied elev rasters
  ##          handle no feature file found
  
  ####    set up cropping mask file and vector of states/etc
  mapcrop <- mapMask(USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                     USParkvec=USParkvec,worldCountryvec=worldCountryvec,
                     mapWindow=mapWindow,
                     mapbuffer=mapbuffer,mapmergebuffer=mapmergebuffer,
                     parkdir=parkDir,
                     workProj4=workProj4)
  mapRectangle <- !is.null(mapWindow)
  statesInMap <- union(expandRegions(unique(toupper(USStatevec)),"US"),
                       expandRegions(unique(toupper(CAProvincevec)),"CANADA")
  )
  if (is.null(statesInMap) & elevDataSource=="Raster") 
    stop(paste0("no states or provinces specified for loading elevations"))
  if (is.null(statesInMap) & featureDataSource=="Raster") 
    stop(paste0("no states or provinces specified for loading features"))
  
  #   now crop to the cropbox 
  if (!is.null(cropbox)) {
    if (writeElevFile | writeFeatureFile | writeShapefiles) 
      warning("cropping map when saving raster/shapefiles.")
    CP <- as(cropbox, "SpatialPolygons")
    sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))
    mapcrop <- rgeos::gIntersection(mapcrop, CP, byid=TRUE)
  }
  plot(mapcrop)  #  which has CRS = workProj4
  
  if (elevDataSource=="Raster") {
    if (!is.null(rasterFileSetNames)) {
      savedNameVec <- rasterFileSetNames
    } else {
      savedNameVec <- statesInMap
    }
    elevations <- loadSavedElevData(savedNameVec=savedNameVec,
                                    rasterDir=rasterDir) 
    if (!is.null(mapWindow))
        elevations <- quickmask(elevations,mapcrop,rectangle=mapRectangle)
  }  else {
    elevations <- loadMapElevData(mapcrop=mapcrop,
                                  mapDataDir=mapDataDir,resstr=resstr)
    #  data loaded from SRTM files is masked
  }
  
  if (writeElevFile) {
    if (!is.null(rasterFileSetWriteNames)) {
      fname <- rasterFileSetWriteNames
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    writeElevRaster(elevations,maxrastercells,rasterDir,
                                     fname=fname)
  }
  ##  Elevation data set up, now load feature raster or at least get shapefiles ready
  ##
  
  featureStack <- NULL
  if (featureDataSource=="Raster") {
    if (!is.null(rasterFileSetNames)) {
      savedNameVec <- rasterFileSetNames
    } else {
      savedNameVec <- statesInMap
    }
    featureStack <- loadSavedFeatureData(savedNameVec=savedNameVec,
                                         rasterDir=rasterDir) %>%
      raster::crop(.,mapcrop) 
  } else if (featureDataSource %in% c("Shapefiles","TIGER")) {
    tmp <- loadShapeFiles(USStatevec,CAProvincevec,mapcrop,
                          shapefileDir,writeShapefiles,
                          shapefileSource=featureDataSource,
                          includeAllRoads=includeAllRoads)
    spTown <- tmp[["spTown"]]
    spRoads <- tmp[["spRoads"]]
    spWaterA <- tmp[["spWaterA"]]
    spWaterL <- tmp[["spWaterL"]]
    #if (!is.null(spTown)) plot(spTown)
    #if (!is.null(spRoads)) plot(spRoads)
    #if (!is.null(spWaterA)) plot(spWaterA)
    #if (!is.null(spWaterL)) plot(spWaterL)
  } else {
    spTown <- NULL
    spRoads <- NULL
    spWaterA <- NULL
    spWaterL <- NULL
  }
  if (writeFeatureFile) { 
    if (!is.null(rasterFileSetWriteNames)) {
      fname <- rasterFileSetWriteNames
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    print("writing feature raster")
    writeFeatureRaster(featureStack,elevations,mapshape=mapcrop,
                       maxrastercells=maxrastercells,
                       rasterDir=rasterDir,
                       fname=fname,
                       spTown=spTown,spRoads=spRoads,
                       spWaterA=spWaterA,spWaterL=spWaterL,
                       polySimplify=0.0,polyMethod="vis", 
                       polyWeighting=0.85,polySnapInt=0.0001)
    #  only need do this if more than one chunk
    print("loading feature raster")
    featureStack <- loadSavedFeatureData(savedNameVec=fname,
                                         rasterDir=rasterDir) 
  }
  
  if (!is.null(featureStack)) { 
    elevations <- raster::addLayer(elevations,featureStack)
    names(elevations[[1]]) <- "elevations"
  }

  # raster(s) already masked by now

  print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
  sfact <- ceiling(sqrt(elevations@ncols*elevations@nrows/(res3dplot*res3dplot)))
  print(paste0("scaling raster down by a factor of ",sfact))
  if (sfact > 1)
    print(system.time(
      elevations <- raster::aggregate(elevations,fact=sfact,fun=max,
                                      expand=TRUE,na.rm=TRUE)
    )[[3]])
  print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
  
  if (featureDataSource %in% c("Shapefiles","TIGER") &
      !writeFeatureFile){
    #  rasterize at the rendered resolution
    featureStack <- buildFeatureStack(elevations,mapshape=mapcrop,
                                      spTown,spWaterA,spWaterL,spRoads,
                                      polySimplify=0.0,polyMethod="vis", 
                                      polyWeighting=0.85,polySnapInt=0.0001)
    
    elevations <- raster::addLayer(elevations,featureStack)
    names(elevations[[1]]) <- "elevations"
  }  
  plot(elevations[["road"]])
  plot(elevations[["waterA"]])
  if (drawRGL)
    draw3DMapTrack(mapRaster=elevations,
                   featureLevels=featureFilter, #towns,roads,waterA,waterL
                   maxElev=maxElev,vScale=vScale,
                   colors=rglColorScheme,
                   citycolor=citycolor,roadcolor=roadcolor,
                   watercolor=watercolor,glaciercolor=glaciercolor,
                   rglNAcolor=rglNAcolor,rglNegcolor=rglNegcolor,
                   saveRGL=saveRGL,mapoutputdir=mapoutputdir,outputName=outputName) 
  return(NULL)
}
