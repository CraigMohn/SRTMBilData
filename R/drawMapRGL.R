drawMapRGL <- function(mapWindow=NULL,
                       USStatevec=NULL,CAProvincevec=NULL,USParkvec=NULL,
                       worldCountryvec=NULL,
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
                       elevDataSource="SRTM",
                       rasterFileSetNames=NULL,rasterFileSetWriteName=NULL,
                       featureDataSource="Shapefiles",
                       writeElevFile=FALSE,writeFeatureFile=FALSE,
                       writeShapefiles=TRUE,includeAllRoads=FALSE,year=2017,
                       rasterDir=NULL,mapDataDir=NULL,shapefileDir=NULL,parkdir=NULL,
                       resstr="_1arc_v3_bil",
                       mapbuffer=0,mapmergebuffer=0,
                       maxrastercells=250000000,maxRasterize=500000,
                       polySimplify=0.0,polyMethod="vis", 
                       polyWeighting=0.85,polySnapInt=0.0001) {

  #  need to do read rasters write smaller rasters
  subsetRasters <- elevDataSource=="Raster" & 
                   featureDataSource=="Raster" &
                   !is.null(rasterFileSetWriteName)
  
  if (elevDataSource=="Raster" & writeElevFile &
           is.null(rasterFileSetWriteName)) {
    warning("will not overwrite elev rasterfileset that was source")
    writeElevFile <- FALSE
  }
  if (featureDataSource=="Raster" & writeFeatureFile &
           is.null(rasterFileSetWriteName)){
    warning("will not overwrite feature rasterfileset that was source")
    writeFeatureFile <- FALSE
  }

  featureFilter <- c(townLevel,roadLevel,waterALevel,waterLLevel)
  
  #############################################################################
  ####    set up masking/cropping mapshape
  mapshape <- mapMask(USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                     USParkvec=USParkvec,worldCountryvec=worldCountryvec,
                     mapWindow=mapWindow,
                     mapbuffer=mapbuffer,mapmergebuffer=mapmergebuffer,
                     parkdir=parkDir,
                     workProj4=workProj4,
                     year=year)
  mapRectangle <- !is.null(mapWindow)
  statesInMap <- union(expandRegions(unique(toupper(USStatevec)),"US"),
                       expandRegions(unique(toupper(CAProvincevec)),"CANADA")
  )
  if (is.null(statesInMap) & 
      is.null(rasterFileSetNames) & 
      elevDataSource=="Raster") 
    stop(paste0("no state/province or dataset names specified for loading elevations"))
  if (is.null(statesInMap) & 
      is.null(rasterFileSetNames) & 
      featureDataSource=="Raster") 
    stop(paste0("no state/province or dataset names specified for loading features"))
  
  #   now crop to the cropbox 
  if (!is.null(cropbox)) {
    if (writeElevFile | writeFeatureFile | writeShapefiles) 
      warning("cropping map when saving raster/shapefiles.")
    CP <- as(cropbox, "SpatialPolygons")
    sp::proj4string(CP) <- CRS(sp::proj4string(mapshape))
    mapshape <- rgeos::gIntersection(mapshape, CP, byid=TRUE)
  }
  plot(mapshape)  #  which has CRS = workProj4
  
  #############################################################################
  ###    load or build elevations raster masked to mapshape
  if (elevDataSource=="Raster") {
    if (!is.null(rasterFileSetNames)) {
      savedNameVec <- rasterFileSetNames
    } else {
      savedNameVec <- statesInMap
    }
    elevations <- loadSavedElevData(savedNameVec=savedNameVec,
                                    rasterDir=rasterDir) 
    if (!is.null(mapWindow))
        elevations <- quickmask(elevations,mapshape,rectangle=mapRectangle)
  }  else {
    elevations <- loadMapElevData(mapshape=mapshape,
                                  mapDataDir=mapDataDir,resstr=resstr)
    #  data loaded from SRTM files is masked
  }
  ###   and write out the elevation raster if requested
  if (writeElevFile) {
    if (!is.null(rasterFileSetWriteName)) {
      fname <- rasterFileSetWriteName
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    writeElevRaster(elevations,maxrastercells,rasterDir,
                                     fname=fname)
    if (writeFeatureFile & !subsetRasters)  
      featuresForElevations(rasterFileSetName=fname,
                            rasterDir=rasterDir,shapefileDir=shapefileDir,
                            USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                            featureDataSource=featureDataSource,
                            writeShapefiles=writeShapefiles,
                            includeAllRoads=includeAllRoads,
                            workProj4=workProj4,
                            maxRasterize=maxRasterize,
                            polySimplify=polySimplify,polyMethod=polymethod, 
                            polyWeighting=polyWeighting,
                            polySnapInt=polySnapInt)
  }
  
  featureStack <- NULL
  if ((featureDataSource=="Raster") | 
      (writeFeatureFile & writeElevFile & !subsetRasters)){
    if (!is.null(rasterFileSetNames)) {
      savedNameVec <- rasterFileSetNames
    } else {
      savedNameVec <- statesInMap
    }
    featureStack <- loadSavedFeatureData(savedNameVec=savedNameVec,
                                         rasterDir=rasterDir) %>%
      raster::crop(.,mapshape) 
    if (subsetRasters) 
      writeFeatureRaster(featureStack,maxrastercells,rasterDir,fname=fname)
  } else if (elevDataSource=="Raster" & writeFeatureFile) {
    if (!is.null(rasterFileSetWriteName)) {
      fname <- rasterFileSetWriteName
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    featuresForElevations(rasterFileSetName=fname,
                          rasterDir=rasterDir,shapefileDir=shapefileDir,
                          USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                          featureDataSource=featureDataSource,
                          writeShapefiles=writeShapefiles,
                          includeAllRoads=includeAllRoads,
                          workProj4=workProj4,
                          maxRasterize=maxRasterize,
                          polySimplify=polySimplify,polyMethod=polymethod, 
                          polyWeighting=polyWeighting,
                          polySnapInt=polySnapInt)
    featureStack <- loadSavedFeatureData(savedNameVec=savedNameVec,
                                       rasterDir=rasterDir)
  } else if (featureDataSource %in% c("Shapefiles","TIGER")) {
    tmp <- loadShapeFiles(USStatevec,CAProvincevec,mapshape,
                          shapefileDir,writeShapefiles,
                          shapefileSource=featureDataSource,
                          includeAllRoads=includeAllRoads,
                          year=year)
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
      is.null(featureStack)) {
    #  rasterize at the rendered resolution
    featureStack <- buildFeatureStack(elevations,mapshape=mapshape,
                                      spTown,spWaterA,spWaterL,spRoads,
                                      maxRasterize=maxRasterize,
                                      polySimplify=0.0,polyMethod="vis", 
                                      polyWeighting=0.85,polySnapInt=0.0001)
    
    elevations <- raster::addLayer(elevations,featureStack)
    names(elevations[[1]]) <- "elevations"
  }  
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
