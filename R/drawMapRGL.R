#' render a 3-D map from elevation and data and 
#'    add a set of paths to the plot
#'
#' \code{drawMapRGL} draw 3D map from elevation and feature data,
#'    and add a set of paths to the plot
#'
#' Create a map displaying the tracks specified.  Tracks may be drawn in one
#'    specified color, with each separate track drawn in a different color from
#'    a specified palette, or as a series of closely spaced dots (which may
#'    be prodded into appearing as a line through judicious choice of point
#'    size and alpha) with the color varying based on speed and the chosen
#'    palette
#'
#' @param paths a set of segments which are multilinestrings. a 
#'    data frame or tibble containing at least: position_lat.dd,
#'    position_lon.dd,(or lat,lon)(both numeric,decimal degrees),
#'    segment(numeric)
#' @param mapwindow a vector of 4 numbers which describe the region drawn.
#'    The format is c(lon_min, lon_max, lat_min, lat_max)
#'    
#' @param USStatevec vector of standard 2-letter abbreviations of US states, or
#'    pre-defined regional aggregations of states for use in defining map or
#'    finding features to include in the map
#' @param CAProvincevec vector of standard 2-letter abbreviations of Canadian
#'    provinvces states, or pre-defined regional aggregations of states for 
#'    use in defining map or finding features to include in the map
#' @param USParkvec vector of US National Park names
#' @param parkdir location of downloaded US Park boundary shapefiles 
#'    downloaded from, for example, https://irma.nps.gov/DataStore/ 
#' @param worldCountryvec vector of 3-letter ISO country abbreviations 
#' @param mapbuffer numeric value to expand defined map
#' @param mapmergebuffer numeric value to expand components of the 
#'    map defined before merging, to eliminate gaps from minor
#'    boundary inconsistencies
#' @param cropbox vector of 4 numbers for cropping the map defined above
#'    The format is c(lon_min, lon_max, lat_min, lat_max).
#' @param rectangularMapdraw draw the entire rectangle enclosing the specified
#'    map area, as long as all relevant state/province rasters have been loaded
#'
#' @param elevDataSource "SRTM" to load data from saved SRTM data, 
#'    "Raster" to load saved raster files
#' @param mapDataDir directory where zipped SRTM data files reside
#' @param resstr suffix on SRTM data files after lon/lat info
#' @param rasterFileSetNames vector of names of saved raster data files
#' 
#' @param featureDataSource "Shapefiles" to load saved shapefiles, "TIGER"
#'    to fetch TIGER data for US states, "Raster" to load saved raster data
#'    from directory specified 
#' @param shapefileDir location to load/save shapefiles
#' @param writeShapefiles write/overwrite shapefiles if TIGER data is used
#' @param year numeric year to use in calls for map boundaries/features
#' @param includeAllRoads include all roads in shapefile, not just highways
#' @param zeroBufferTowns use zero buffer trick to repair town polygon shapefile
#' @param zeroBufferWater use zero buffer trick to repair water polygon shapefile

#' @param writeElevFile save the elevation raster files 
#' @param writeFeatureFile save the feature data raster files
#' @param rasterDir location to load and save raster files
#' @param rasterFileSetWriteName name to use for the raster fileset written

#' @param drawRGL if TRUE, open a window with the interactive 3D map
#' @param res3dplot used to scale resolution, the maximum size the square root
#'    of the output height times width.  aggregation done in integer multiples.
#' @param maxElev all elevations greater than this are colored the same
#' @param vScale vertical scale parameter, use larger for smaller areas
#' @param townLevel display towns ranked this number or higher: 
#'    3=all towns     5=larger towns (in US >50k)
#' @param roadLevel display roads ranked this number or higher:
#'    2=Service Drive, Bike Path, etc      3=Local Street
#'    4=Secondary Hwy                      5=Primary Hwy/Transit
#' @param waterALevel display areal water ranked this number or higher:
#'    2=res/treatmentpond/pit/quarry       3=lake/pond/swamp/stream
#'    4=class 2 or 3 bigger than 1k ha     5=named lake/pond/swamp/stream
#'    6=large lake/pond/swamp/stream       7=Ocean/Bay/Est/Sound
#'    8=glacier
#' @param waterLLevel display linear water ranked this number or higher:
#'    2=canal/ditch                        3=braided stream
#'    4=stream/river                       5=named stream/river
#'    6=named stream/river containing the string "RIV"
#' @param rglColorScheme name of color scheme from 
#'     c("default","beach","viridis","plasma","terrain","oleron","snow","oslo",
#'       "desert","lajolla","niccoli","bright",
#'       "bing","maptoolkit-topo","nps","apple-iphoto")
#' @param mapColorDepth number of bits per color channel in map image
#' @param rglNAcolor color used to display NA elevations
#' @param rglNegcolor color used to display elevations below sea level
#' @param citycolor color used to display cities
#' @param watercolor color used to display water,
#' @param roadcolor color used to display roads
#' @param glaciercolor color used to display glaciers
#' @param rglShininess number controlling surface shininess in rgl rendering
#' @param rglSpecular light color for specular light
#' @param rglDiffuse light color for diffuse light
#' @param rglAmbient light color for ambient light
#' @param rglSmooth use Gouraud shading if TRUE
#' @param rglAlpha alpha value for surface
#' @param rglAntiAlias antialias points and lines when drawing surface
#' @param rglTheta coordinate for light source
#' @param rglPhi coordinate for light source
#' @param saveRGL save the map to an html file
#' @param mapoutputdir location for saved html map file
#' @param outputName name of saved html map

#' @param workProj4 coordinte reference projection string
#' @param maxrastercells maximum number of cells in each written raster
#' @param maxRasterize number of items for calls to velox$rasterize
#' @param polySimplify amount of polygon simplification, see help
#'    for rmapshaper::ms_simplify 
#' @param polyMethod simplification method either "vis" or "dp"
#' @param polyWeighting see help for rmapshaper::ms_simplify 
#' @param polySnapInt see help for rmapshaper::ms_simplify
  

#' @return NULL
#'
#' @export
drawMapRGL <- function(paths=NULL,
                       # map shape parameters
                       mapWindow=NULL,
                       USStatevec=NULL,CAProvincevec=NULL,
                       USParkvec=NULL,parkdir=NULL,
                       worldCountryvec=NULL,
                       cropbox=NULL,mapbuffer=0,mapmergebuffer=0,
                       rectangularMap=TRUE,
                       #  elevation data parameters
                       elevDataSource="SRTM",
                       mapDataDir=NULL,resstr="_1arc_v3_bil",
                       rasterFileSetNames=NULL,
                       #  water/town/road features
                       featureDataSource="Shapefiles",
                       shapefileDir=NULL,writeShapefiles=TRUE,
                       year=2017,includeAllRoads=FALSE,
                       zeroBufferTowns=FALSE,
                       zeroBufferWater=FALSE,
                       #  raster save and location parameters
                       writeElevFile=FALSE,writeFeatureFile=FALSE,
                       rasterDir=NULL,rasterFileSetWriteName=NULL,
                       #  plotting control
                       drawRGL=TRUE,
                       res3dplot=2500,maxElev=3000,vScale=1.5,
                       townLevel=3,roadLevel=4,waterALevel=4,waterLLevel=5,
                       rglColorScheme="default",mapColorDepth=16,
                       rglNAcolor="Blue",rglNegcolor="Red",
                       citycolor="SlateGray",watercolor="Blue",
                       roadcolor="Black",glaciercolor="White",
                       rglShininess=0,
                       rglSpecular="black", rglDiffuse="white", 
                       rglAmbient="white", rglEmission="black",
                       rglSmooth=TRUE,
                       rglAlpha=1.0,
                       rglAntiAlias=TRUE,
                       rglTheta=0,
                       rglPhi=15,
                       saveRGL=FALSE,mapoutputdir=NULL,outputName=NULL,
                       #  CRS, rasterization control
                       workProj4="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 +no_defs",
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
  mapshape <- mapMask(USStatevec=USStatevec,
                      CAProvincevec=CAProvincevec,
                      USParkvec=USParkvec,
                      worldCountryvec=worldCountryvec,
                      mapWindow=mapWindow,
                      mapbuffer=mapbuffer,
                      mapmergebuffer=mapmergebuffer,
                      parkDir=parkDir,workProj4=workProj4,
                      year=year)
  statesInMap <- union(expandRegions(unique(toupper(USStatevec)),"US"),
                       expandRegions(unique(toupper(CAProvincevec)),"CANADA"))
  
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
    CP <- as(raster::extent(cropbox), "SpatialPolygons")
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
    elevations <- loadSavedElevData(savedNameVec,rasterDir) 
    if (!is.null(mapWindow) | !is.null(USParkvec) | !is.null(cropbox)) {
      print("cropping")
      print(paste0("  ",round(system.time(
        elevations <- quickmask(elevations,mapshape,rectangle=TRUE)
      )[[3]],digits=2)))
    }
  }  else {
    elevations <- loadMapElevData(mapshape,mapDataDir,resstr)
    #  data loaded from SRTM files is masked
  }
  ###   and write out the elevation raster if requested
  if (writeElevFile) {
    if (!is.null(rasterFileSetWriteName)) {
      fname <- rasterFileSetWriteName
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    writeElevRaster(elevations=elevations,
                    maxrastercells=maxrastercells,
                    rasterDir=rasterDir,
                    fname=fname)
    if (writeFeatureFile & !subsetRasters)  
      featuresForElevations(rasterFileSetName=fname,
                            rasterDir=rasterDir,
                            shapefileDir=shapefileDir,
                            USStatevec=USStatevec,
                            CAProvincevec=CAProvincevec,
                            featureDataSource=featureDataSource,
                            writeShapefiles=writeShapefiles,
                            includeAllRoads=includeAllROads,
                            zeroBufferWater=zeroBufferWater,
                            workProj4=workProj4,
                            maxRasterize=maxRasterize,
                            polySimplify=polySimplify,
                            polyMethod=polyMethod, 
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
    featureStack <- loadSavedFeatureData(savedNameVec,rasterDir)
    if (!is.null(mapWindow) | !is.null(USParkvec) | !is.null(cropbox)) {
      print("cropping")
      print(paste0("  ",round(system.time(
        featureStack <- quickmask(featureStack,mapshape,rectangle=TRUE)
      )[[3]],digits=2)))
    }
    if (subsetRasters) 
      writeFeatureRaster(featureStack=featureStack,
                         maxrastercells=maxrastercells,
                         rasterDir=rasterDir,
                         fname=fname)
  } else if (elevDataSource=="Raster" & writeFeatureFile) {
    if (!is.null(rasterFileSetWriteName)) {
      fname <- rasterFileSetWriteName
    } else {
      fname <- paste0(statesInMap,collapse="")
    }
    featuresForElevations(rasterFileSetName=fname,
                          rasterDir=rasterDir,
                          shapefileDir=shapefileDir,
                          USStatevec=USStatevec,
                          CAProvincevec=CAProvincevec,
                          featureDataSource=featureDataSource,
                          writeShapefiles=writeShapefiles,
                          includeAllRoads=includeAllROads,
                          year=year,
                          zeroBufferTowns=zeroBufferTowns,
                          zeroBufferWater=zeroBufferWater,
                          workProj4=workProj4,
                          maxRasterize=maxRasterize,
                          polySimplify=polySimplify,
                          polyMethod=polyMethod, 
                          polyWeighting=polyWeighting,
                          polySnapInt=polySnapInt)
    featureStack <- loadSavedFeatureData(savedNameVec,rasterDir)
  } else if (featureDataSource %in% c("Shapefiles","TIGER")) {
    tmp <- loadShapeFiles(USStatevec=USStatevec,
                          CAProvincevec=CAProvincevec,
                          mapshape=mapshape,
                          shapefileDir=shapefileDir,
                          writeShapefiles=writeShapefiles,
                          shapefileSource=featureDataSource,
                          includeAllRoads=includeAllROads,
                          year=year,
                          zeroBufferTowns=zeroBufferTowns,
                          zeroBufferWater=zeroBufferWater)
    spTown <- tmp[["spTown"]]
    spRoads <- tmp[["spRoads"]]
    spWaterA <- tmp[["spWaterA"]]
    spWaterL <- tmp[["spWaterL"]]
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
  ##  rasters are cropped but not necessarily masked by this point
  print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
  sfact <- ceiling(sqrt((as.double(elevations@ncols)/res3dplot)*
                        (as.double(elevations@nrows)/res3dplot)))
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
    featureStack <- buildFeatureStack(baseLayer=elevations,
                                      mapshape=mapshape,
                                      spTown=spTown,
                                      spWaterA=spWaterA,
                                      spWaterL=spWaterL,
                                      spRoads=spRoads,
                                      filterVec=featureFilter,
                                      maxRasterize=maxRasterize,
                                      polySimplify=polySimplify,
                                      polyMethod=polyMethod, 
                                      polyWeighting=polyWeighting,
                                      polySnapInt=polySnapInt)
    
    elevations <- raster::addLayer(elevations,featureStack)
    names(elevations[[1]]) <- "elevations"
  }  
  ##  need to mask if not a rectangle - whole states/provinces/countries are masked already,
  ##     so only need to worry about Parks
  if (is.null(mapWindow) & !rectangularMap & !is.null(USParkvec))
    print(paste0("  ",round(system.time(
      elevations <- quickmask(elevations,mapshape,rectangle=FALSE)
    )[[3]],digits=2)," rasterizing"))

  if (drawRGL)
    draw3DMapTrack(mapRaster=elevations,trackdf=paths,
                   featureLevels=featureFilter, #towns,roads,waterA,waterL
                   maxElev=maxElev,vScale=vScale,
                   rglColorScheme=rglColorScheme,
                   mapColorDepth=mapColorDepth,
                   citycolor=citycolor,roadcolor=roadcolor,
                   watercolor=watercolor,glaciercolor=glaciercolor,
                   rglNAcolor=rglNAcolor,rglNegcolor=rglNegcolor,
                   rglShininess=rglShininess,rglSmooth=rglSmooth,
                   rglAlpha=rglAlpha,rglAntiAlias=rglAntiAlias,
                   rglSpecular=rglSpecular,rglDiffuse=rglDiffuse,
                   rglAmbient=rglAmbient,rglEmission=rglEmission,
                   rglTheta=rglTheta,rglPhi=rglPhi,saveRGL=saveRGL,
                   mapoutputdir=mapoutputdir,outputName=outputName) 
  return(NULL)
}
