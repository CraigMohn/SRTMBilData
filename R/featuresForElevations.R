elevationsToRaster <- function(rasterFileSetName="default",
                               USStatevec=NULL,CAProvincevec=NULL,
                               USParkvec=NULL,worldCountryvec=NULL,
                               mapWindow=NULL,cropbox=NULL,
                               rasterDir=NULL,mapDataDir=NULL,parkdir=NULL,
                               mapbuffer=0,mapmergebuffer=0,
                               maxrastercells=250000000,
                               workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                               resstr="_1arc_v3_bil") {
  
  mapshape <- mapMask(USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                      USParkvec=USParkvec,worldCountryvec=worldCountryvec,
                      mapWindow=mapWindow,
                      mapbuffer=mapbuffer,mapmergebuffer=mapmergebuffer,
                      parkdir=parkDir,
                      workProj4=workProj4)
  
  #   now crop to the cropbox 
  if (!is.null(cropbox)) {
    warning("cropping map when saving raster/shapefiles.")
    CP <- as(cropbox, "SpatialPolygons")
    sp::proj4string(CP) <- CRS(sp::proj4string(mapshape))
    mapshape <- rgeos::gIntersection(mapshape, CP, byid=TRUE)
  }
  plot(mapshape)  #  which has CRS = workProj4
  
  elevations <- loadMapElevData(mapshape=mapshape,
                                mapDataDir=mapDataDir,resstr=resstr)
  writeElevRaster(elevations,maxrastercells,rasterDir,
                  fname=rasterFileSetName)
  return(elevations)
}
featuresForElevations <- function(rasterFileSetName,
                                  rasterDir=NULL,shapefileDir=NULL,
                                  USStatevec=NULL,
                                  CAProvincevec=NULL,
                                  featureDataSource="Shapefiles",
                                  writeShapefiles=TRUE,includeAllRoads=FALSE,
                                  workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                                  mapbuffer=0,mapmergebuffer=0,
                                  maxRasterize=100000,
                                  polySimplify=0.0,polyMethod="vis", 
                                  polyWeighting=0.85,polySnapInt=0.0001) {

  #   build mapshape
  mapshape <- mapMask(USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                      USParkvec=NULL,worldCountryvec=NULL,
                      mapWindow=NULL,
                      mapbuffer=mapbuffer,mapmergebuffer=mapmergebuffer,
                      parkdir=NULL,
                      workProj4=workProj4)
  plot(mapshape)

  #   get shapefiles from State and Province vecs
  tmp <- loadShapeFiles(USStatevec,CAProvincevec,mapshape=mapshape,
                        shapefileDir=shapefileDir,
                        writeShapefiles=writeShapefiles,
                        shapefileSource=featureDataSource,
                        includeAllRoads=includeAllRoads)
  spTown <- tmp[["spTown"]]
  spRoads <- tmp[["spRoads"]]
  spWaterA <- tmp[["spWaterA"]]
  spWaterL <- tmp[["spWaterL"]]
  

  fvec <- list.files(path=paste0(rasterDir,"/",rasterFileSetName),
                     pattern=paste0(rasterFileSetName,"elevs[0-9]{,2}.grd"))
  for (fn in fvec) {
    fname <- sub(paste0(rasterFileSetName,"elevs"),
                 paste0(rasterFileSetName,"features"),
                 sub(".grd","",fn))
    print(paste0("loading ",fn))
    elevations <- raster(paste0(rasterDir,"/",rasterFileSetName,"/",fn))
    featureStack <- buildFeatureStack(elevations,mapshape=mapshape,
                                      spTown=spTown,spWaterA=spWaterA,
                                      spWaterL=spWaterL,spRoads=spRoads,
                                      maxRasterize=maxRasterize,
                                      polySimplify=polySimplify,
                                      polyMethod=polyMethod, 
                                      polyWeighting=polyWeighting,
                                      polySnapInt=polySnapInt)
    print(paste0("writing ",fname,".grd"))
    writeRaster(featureStack,
                file=paste0(rasterDir,"/",rasterFileSetName,"/",
                            fname,
                            ".grd"),
                bylayer=TRUE,suffix="names",   
                datatype="INT1S",overwrite=TRUE)   
    featureStack <- NULL
    gc()
  }
  return(NULL)
}
  