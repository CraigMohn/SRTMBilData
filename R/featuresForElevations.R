#' Create and save a raster fileset specified 
#'
#' \code{elevationsToRaster} create a raster fileset of elevations
#'
#' Create a raster file set containing elevation data from SRTM BIL 
#'   zip files stored in a local directory
#'
#' @param rasterFileSetName names of saved raster data files
#' @param rasterDir location to load and save raster files
#' @param mapDataDir directory where zipped SRTM data files reside
#' @param USStatevec vector of standard 2-letter abbreviations of US states, or
#'    pre-defined regional aggregations of states for use in defining map or
#'    finding features to include in the map
#' @param CAProvincevec vector of standard 2-letter abbreviations of Canadian
#'    provinvces states, or pre-defined regional aggregations of states for 
#'    use in defining map or finding features to include in the map
#' @param USParkvec vector of US National Park names
#' @param worldCountryvec vector of 3-letter ISO country abbreviations 
#' @param mapwindow a vector of 4 numbers which describe the region drawn.
#'    The format is c(lon_min, lon_max, lat_min, lat_max)
#' @param cropbox vector of 4 numbers for cropping the map defined above
#'    The format is c(lon_min, lon_max, lat_min, lat_max).
#' @param parkdir location of downloaded US Park boundary shapefiles 
#'    downloaded from, for example, https://irma.nps.gov/DataStore/ 
#' @param mapbuffer numeric value to expand defined map
#' @param mapmergebuffer numeric value to expand components of the 
#'    map defined before merging, to eliminate gaps from minor
#'    boundary inconsistencies
#' @param maxrastercells maximum number of cells in each written raster
#' @param workProj4 coordinte reference projection string
#' @param year numeric year to use in calls for map boundaries/features
#' @param resstr suffix on SRTM data files after lon/lat info

#' @return a rasterLayer containing the elevations
#'
#' @export
elevationsToRaster <- function(rasterFileSetName="default",
                               rasterDir,mapDataDir,
                               USStatevec=NULL,CAProvincevec=NULL,
                               USParkvec=NULL,worldCountryvec=NULL,
                               mapWindow=NULL,cropbox=NULL,
                               parkdir=NULL,
                               mapbuffer=0,mapmergebuffer=0,
                               maxrastercells=250000000,
                               workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                               year=2017,
                               resstr="_1arc_v3_bil") {
  
  
  return(drawMapRGL(paths=NULL,
                    mapWindow,USStatevec,CAProvincevec,
                    USParkvec,parkdir,
                    worldCountryvec,
                    cropbox,mapbuffer,mapmergebuffer,
                    rectangularMap=FALSE,drawRGL=FALSE,
                    elevDataSource="SRTM",mapDataDir,resstr,
                                rasterFileSetNames=NULL,
                    featureDataSource="none",
                    year,WriteElevFile=TRUE,writeFeatureFile=FALSE,
                    rasterDir,workProj4,
                    maxrastercells))
  
}
#' create and save the laters of a rasterStack containing the feature data for 
#'   an elevation raster
#'
#' \code{featuresForElevations} create and save the laters of a rasterStack 
#'   containing the feature data for an elevation raster
#'
#' @param rasterFileSetNames vector of names of saved raster data files
#' @param rasterDir location to load and save raster files
#' @param shapefileDir location to load/save shapefiles
#' @param USStatevec vector of standard 2-letter abbreviations of US states, or
#'    pre-defined regional aggregations of states for use in defining map or
#'    finding features to include in the map
#' @param CAProvincevec vector of standard 2-letter abbreviations of Canadian
#'    provinvces states, or pre-defined regional aggregations of states for 
#'    use in defining map or finding features to include in the map
#' @param featureDataSource "Shapefiles" to load saved shapefiles, "TIGER"
#'    to fetch TIGER data for US states, "Raster" to load saved raster data
#'    from directory specified 
#' @param writeShapefiles write/overwrite shapefiles if TIGER data is used
#' @param includeAllRoads include all roads in shapefile, not just highways
#' @param year numeric year to use in calls for map boundaries/features
#' @param zeroBufferTowns use zero buffer trick to repair town polygon shapefile
#' @param zeroBufferWater use zero buffer trick to repair water polygon shapefile
#' @param workProj4 coordinte reference projection string
#' @param mapbuffer numeric value to expand defined map
#' @param mapmergebuffer numeric value to expand components of the 
#'    map defined before merging, to eliminate gaps from minor
#'    boundary inconsistencies
#' @param maxRasterize number of items for calls to velox$rasterize
#' @param polySimplify amount of polygon simplification, see help
#'    for rmapshaper::ms_simplify 
#' @param polyMethod simplification method either "vis" or "dp"
#' @param polyWeighting see help for rmapshaper::ms_simplify 
#' @param polySnapInt see help for rmapshaper::ms_simplify
#'
#' @return NULL
#'
#' @export
featuresForElevations <- function(rasterFileSetName,
                                  rasterDir,shapefileDir=NULL,
                                  USStatevec=NULL,
                                  CAProvincevec=NULL,
                                  featureDataSource="Shapefiles",
                                  writeShapefiles=TRUE,includeAllRoads=FALSE,
                                  year=2017,
                                  zeroBufferTowns=FALSE,
                                  zeroBufferWater=FALSE,
                                  workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                                  mapbuffer=0,mapmergebuffer=0,
                                  maxRasterize=100000,
                                  polySimplify=0.0,polyMethod="vis", 
                                  polyWeighting=0.85,polySnapInt=0.0001,...) {
  #   build mapshape
  mapshape <- mapMask(USStatevec=USStatevec,CAProvincevec=CAProvincevec,
                      USParkvec=NULL,worldCountryvec=NULL,
                      mapWindow=NULL,
                      mapbuffer=mapbuffer,mapmergebuffer=mapmergebuffer,
                      parkDir=NULL,
                      workProj4=workProj4,year=year)
  plot(mapshape)
  #   get shapefiles from State and Province vecs
  tmp <- loadShapeFiles(USStatevec,CAProvincevec,mapshape=mapshape,
                        shapefileDir=shapefileDir,
                        writeShapefiles=writeShapefiles,
                        shapefileSource=featureDataSource,
                        includeAllRoads=includeAllRoads,
                        year=year,
                        zeroBufferTowns=zeroBufferTowns,
                        zeroBufferWater=zeroBufferWater)
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
  