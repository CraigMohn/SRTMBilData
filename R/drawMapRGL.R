drawMapRGL <- function(USStatevec=NULL,CAProvincevec=NULL,USParkvec=NULL,
                       worldCountryvec=NULL,mapWindow=NULL,
                       routeSL=NULL,cropbox=NULL,
                       res3dplot=3000,maxElev=3000,vScale=1.5,
                       townLevel=3,roadLevel=3,waterALevel=4,waterLLevel=5,
                       rglColorScheme="default",
                       rglNAcolor="Blue",rglNegcolor="Red",
                       citycolor="Magenta",watercolor="Blue",
                       roadcolor="Black",glaciercolor="White",
                       drawRGL=TRUE,
                       saveRGL=FALSE,mapoutputdir=NULL,outputName=NULL,
                       workProj4="+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                       elevDataSource="Raster",useRasterFileSets=NULL,
                       featureDataSource="Shapefiles",
                       writeElevFile=FALSE,writeFeatureFile=FALSE,
                       writeShapefiles=FALSE,
                       rasterDir=NULL,mapDataDir=NULL,shapefileDir=NULL,parkdir=NULL,
                       resstr="_1arc_v3_bil",
                       mapbuffer=0,mapmergebuffer=0,
                       maxrastercells=500000000,
                       polySimplify=0.0,polyMethod="vis", 
                       polyWeighting=0.85,polySnapInt=0.0001) {

  if (elevDataSource=="Raster" & writeElevFile)
    stop("don't ask me to read the raster elevation data and then write it")
  if (featureDataSource=="Raster" & writeFeatureFile)
    stop("don't ask me to read the raster feature data and then write it")
  if (featureDataSource=="Shapefiles" & writeShapefiles)
    stop("don't ask me to read the shapefile data and then write it")
  
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
    elevations <- loadSavedElevData(savedNameVec=statesInMap,
                                    rasterDir=rasterDir) 
  }  else {
    elevations <- loadMapElevData(mapcrop=mapcrop,
                                  mapDataDir=mapDataDir,resstr=resstr)
  }
  if (writeElevFile) writeElevRaster(elevations,maxrastercells,rasterDir,
                                     fname=paste0(statesInMap,collapse=""))
  
  if (featureDataSource=="Raster") {
    featureStack <- loadSavedFeatureData(savedNameVec=statesInMap,
                                         rasterDir=rasterDir) 
  } else {
    tmp <- loadShapeFiles(USStatevec,CAProvincevec,mapcrop,
                          shapefileDir,writeShapefiles,
                          shapefileSource=featureDataSource)
    spTown <- tmp[["spTown"]]
    spRoads <- tmp[["spRoads"]]
    spWaterA <- tmp[["spWaterA"]]
    spWaterL <- tmp[["spWaterL"]]
    featureStack <- buildFeatureStack(elevations,mapshape=mapcrop,
                                      spTown,spWaterA,spWaterL,spRoads,
                                      polySimplify=0.0,polyMethod="vis", 
                                      polyWeighting=0.85,polySnapInt=0.0001)
  }
  if (writeFeatureFile) 
    writeFeatureRaster(featureStack,
                       maxrastercells=maxrastercells,
                       rasterDir=rasterDir,
                       fname=paste0(statesInMap,collapse=""))
  
  elevations <- raster::addLayer(elevations,featureStack)
  names(elevations[[1]]) <- "elevations"
  
  # crop raster now since it's relatively fast
  elevations <- raster::crop(elevations,mapcrop)
  
  print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
  sfact <- ceiling(sqrt(elevations@ncols*elevations@nrows/(res3dplot*res3dplot)))
  print(paste0("scaling raster down by a factor of ",sfact))
  if (sfact > 1)
    print(system.time(
      elevations <- raster::aggregate(elevations,fact=sfact,fun=max,
                                      expand=TRUE,na.rm=TRUE)
    )[[3]])
  print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
  
  # mask raster after aggregation for speed
  elevations <- raster::mask(elevations,mapcrop)
  
  if (drawRGL)
    draw3DMapRoute(mapRaster=elevations,
                   featureLevels=featureFilter, #towns,roads,waterA,waterL
                   maxElev=maxElev,vScale=vScale,
                   colors=rglColorScheme,
                   citycolor=citycolor,roadcolor=roadcolor,
                   watercolor=watercolor,glaciercolor=glaciercolor,
                   rglNAcolor=rglNAcolor,rglNegcolor=rglNegcolor,
                   saveRGL=saveRGL,mapoutputdir=mapoutputdir,outputName=outputName) 
  return(NULL)
}
