loadShapeFiles <- function(USStatevec,CAProvincevec,mapshape,
                           shapefileDir,writeShapefiles,
                           shapefileSource="Shapefiles",
                           includeAllRoads=FALSE,
                           zeroBufferTowns=FALSE,zeroBufferWater=FALSE,
                           year=2017) {
  workCRS <- raster::crs(mapshape)
  spList <- list()
  if (!is.null(USStatevec)) {
    tmp <- USFeatures(USStatevec=expandRegions(unique(toupper(USStatevec)),"US"),
                      workCRS=workCRS,
                      shapefileSource=shapefileSource,
                      shapefileDir=shapefileDir,
                      writeShapefiles=writeShapefiles,
                      includeAllRoads=includeAllRoads,
                      year=year)
    spList <- rbindList(spList,tmp)
  }
  if (!is.null(CAProvincevec)) {
    tmp <- CAFeatures(CAProvincevec=expandRegions(unique(toupper(CAProvincevec)),"CANADA"),
                      workCRS=workCRS,
                      shapefileDir=shapefileDir)
    spList <- rbindList(spList,tmp)
  }
  if (zeroBufferTowns & !is.null(spList[["spTown"]])) {
    spTown <- spList[["spTown"]]
    bufferProj4 <- UTMProj4(lon=(extent(spTown)[1]+extent(spTown)[2])/2,
                            lat=(extent(spTown)[3]+extent(spTown)[4])/2)
    spTown <- rgeos::gBuffer(sp::spTransform(spTown, CRS(bufferProj4)),
                             width=0,byid=TRUE)
    spList[["spTown"]] <- sp::spTransform(spTown , workCRS)
  }
  if (zeroBufferWater & !is.null(spList[["spWaterA"]])) {
    spWaterA <- spList[["spWaterA"]]
    bufferProj4 <- UTMProj4(lon=(extent(spWaterA)[1]+extent(spWaterA)[2])/2,
                            lat=(extent(spWaterA)[3]+extent(spWaterA)[4])/2)
    spWaterA <- rgeos::gBuffer(sp::spTransform(spWaterA, CRS(bufferProj4)),
                               width=0,byid=TRUE)
    spList[["spWaterA"]] <- sp::spTransform(spWaterA , workCRS)
  }
  #  Don't need to check boundaries - 
  #   saved chunks were appropriately masked/filtered and 
  #   TIGER is masked if needed when fetched
  #   canada data was categorized when shapefiles were created
  return(spList)
}
USFeatures <- function(USStatevec,workCRS,
                       shapefileDir,writeShapefiles,
                       shapefileSource="TIGER",
                       includeAllRoads=FALSE,
                       year=2017) {
  if (is.null(USStatevec)) stop("empty USStatevec") 
  if (is.null(shapefileDir) & 
      (!(shapefileSource == "TIGER") | writeShapefiles) )
    stop("need shapefileDir")
  spList <- list()
  for (st in USStatevec) {
    if (shapefileSource == "TIGER" | 
        !file.exists(paste0(shapefileDir,"/",st,"Town.shp"))) {
      tmp <- USTigerFeatures(st,workCRS=workCRS,
                             shapefileDir=shapefileDir,
                             writeShapefiles=writeShapefiles,
                             includeAllRoads=includeAllRoads,
                             year=year)
    } else {
      tmp <- readShapeFiles(st,shapefileDir,workCRS)
    }
    spList <- rbindList(spList,tmp)
  }   
  return(spList)
}
USTigerFeatures <- function(st,workCRS,
                            shapefileDir,writeShapefiles=TRUE,
                            includeAllRoads=FALSE,
                            year=2017) {
  stMask <- sp::spTransform(tigris::counties(st,year=year),workCRS)
  
  spTown <- sp::spTransform(tigris::urban_areas(year=year),workCRS)
  spTown <- spTown[,(names(spTown) %in% c("NAME10","UATYP10"))]
  colnames(spTown@data) <- c("NAME","TYPE")
  spTown <- sxdfMask(spTown,stMask,keepTouch=FALSE) #keep if town touches the state
  # type is   U pop >= 50k, 2.5k <= C pop < 50k   
  
  # spRoads <- sp::spTransform(tigris::primary_secondary_roads(st),workCRS) # SpatialPolygonsDF
  # spRoads <- spRoads[,(names(spRoads) %in% c("FULLNAME","MTFCC"))]
  # colnames(spRoads@data) <- c("NAME","TYPE")  
  # type is   S1100=secondary   S1200=Primary
  
  tmpR <- NULL
  tmpA <- NULL
  tmpL <- NULL
  for (c in tigris::list_counties(st)[["county_code"]]) {
    if (!(c == "515" & st == "VA")) {   #Bedford Town not in data?!?
      #spatialPolygon dataframe
      tmpA <- rbind_NULLok(tmpA,
                           sp::spTransform(tigris::area_water(st,c,year=year),workCRS)) 
      #spatialLines dataframe
      tmpL <- rbind_NULLok(tmpL,
                           sp::spTransform(tigris::linear_water(st,c,year=year),workCRS)) 
      tmpR <- rbind_NULLok(tmpR,
                           sp::spTransform(tigris::roads(st,c,year=year),workCRS)) 
    }
  }  
  tmpA$size <- as.numeric(tmpA$AWATER) + as.numeric(tmpA$ALAND)
  spWaterA <- tmpA[,(names(tmpA) %in% c("FULLNAME","MTFCC","size"))]
  colnames(spWaterA@data) <- c("NAME","TYPE","size")  
  # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
  #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
  
  spWaterL <- tmpL[,(names(tmpL) %in% c("FULLNAME","MTFCC"))]
  colnames(spWaterL@data) <- c("NAME","TYPE")  
  # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
  
  spRoads <- tmpR[,(names(tmpR) %in% c("FULLNAME","MTFCC"))]
  colnames(spRoads@data) <- c("NAME","TYPE")  
  # type is   S1100=Primary        S1200=Secondary      S1400=Local St 
  #           S1500=Vehicular Trl  S1630=Ramp           S1640=Service Drive 
  #           S1710=Walkway        S1720=Stairway       S1730=Alley
  #           S1740=Private Rd     S1750=Internal       S1780=Pkgng Lot 
  #           S1820=Bike Path      S1830=Bridle Trl
  
  if (writeShapefiles) 
    writeShapeFiles(st,shapefileDir=shapefileDir,
                    spTown=spTown,spRoads=spRoads,
                    spWaterA=spWaterA,spWaterL=spWaterL)
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
CAFeatures <- function(CAProvincevec,workCRS,shapefileDir) {
  if (is.null(shapefileDir)) stop("need shapefileDir")
  spList <- list()
  ##  files are already filtered for city/road/polygonwater importance
  for (pr in CAProvincevec) {
    tmp <- readShapeFiles(pr,shapefileDir,workCRS)
    spList <- rbindList(spList,tmp)
  }
  return(spList)
}
readShapeFiles <- function(stname,shapefileDir,workProj4,
                           noisy=FALSE,silent=FALSE) {
  if (!silent) print(paste0("loading shapefiles for ",stname))
  spTown <- sp::spTransform(raster::shapefile(paste0(shapefileDir,"/",stname,
                                                     "Town.shp"))
                            ,workProj4)
  spRoads <- sp::spTransform(raster::shapefile(paste0(shapefileDir,"/",stname,
                                                      "Roads.shp"))
                             ,workProj4)
  spWaterA <- sp::spTransform(raster::shapefile(paste0(shapefileDir,"/",stname,
                                                       "WaterA")),
                              workProj4)
  spWaterL <- sp::spTransform(raster::shapefile(paste0(shapefileDir,"/",stname,
                                                       "WaterL")),
                              workProj4)
  if (!silent) print("done loading shapefiles")
  return(list("spTown"=spTown,"spRoads"=spRoads,
              "spWaterA"=spWaterA,"spWaterL"=spWaterL))
}
