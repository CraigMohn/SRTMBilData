buildFeatureStack <- function(baseLayer,mapshape,
                              spTown,spWaterA,spWaterL,spRoads,
                              filterVec=c(1,1,1,1),  #towns,roads,waterA,waterL
                              maxRasterize=10000,
                              polySimplify=0,polyMethod="vis",
                              polyWeighting=0.85,polySnapInt=0.0001  ) {
  #  baseLayer a rasterLayer
  #  returns a rasterStack  which has 4 layers which can be stored as ints
  bLrect <- as(raster::extent(baseLayer), "SpatialPolygons")
  sp::proj4string(bLrect) <- sp::proj4string(spTown)
  if (!is.null(spTown)) {
    print("towns")
    spTown <- sxdfMask(spTown, bLrect)
  }
  if (!is.null(spTown)) {
    print(paste0(nrow(spTown)," towns to process"))
    spTown$value <- townRank(spTown@data[,"TYPE"],
                             spTown@data[,"NAME"])
    spTown <- spTown[spTown$value>=filterVec[1],]
    tlayer <- shapeToRasterLayer(sxdf=spTown,
                                 templateRaster=baseLayer,
                                 maxRasterize=maxRasterize,
                                 polySimplify=polySimplify,
                                 polyMethod=polyMethod,
                                 polyWeighting=polyWeighting,
                                 polySnapInt=polySnapInt)
    if (!is.finite(tlayer@data@min)) {
      warning("tlayer mess-up")
      print(tlayer@data@min)
    }
  } else {
    tlayer <- raster::raster(nrows=nrow(baseLayer),
                             ncols=ncol(baseLayer),
                             ext=extent(baseLayer),
                             crs=crs(baseLayer),
                             vals=0)
    print("no towns to add")
  }
  if (!is.null(spWaterA)) {
    print("water polygons")
    spWaterA <- sxdfMask(spWaterA,bLrect)
  }
  if (!is.null(spWaterA)) {
    print(paste0(nrow(spWaterA)," water areas to process"))
    spWaterA$value <- waterARank(spWaterA@data[,"TYPE"],
                                  spWaterA@data[,"NAME"],
                                  spWaterA@data[,"size"] )      
    spWaterA <- spWaterA[spWaterA$value>=filterVec[3],]
    wAlayer <- shapeToRasterLayer(sxdf=spWaterA,
                                  templateRaster=baseLayer,
                                  maxRasterize=maxRasterize,
                                  polySimplify=polySimplify,
                                  polyMethod=polyMethod,
                                  polyWeighting=polyWeighting,
                                  polySnapInt=polySnapInt)
    if (!is.finite(wAlayer@data@min)) {
      warning("wAlayer mess-up")
      print(wAlayer@data@min)
    }
  } else {
    wAlayer <- raster::raster(nrows=nrow(baseLayer),
                              ncols=ncol(baseLayer),
                              ext=extent(baseLayer),
                              crs=crs(baseLayer),
                              vals=0)
    print("no water polygons to add")
  }
  if (!is.null(spWaterL)) {
    print("water lines")
    spWaterL <- sxdfMask(spWaterL,bLrect)
  }
  if (!is.null(spWaterL)) {
    print(paste0(nrow(spWaterL)," water lines to process"))
    spWaterL$value <- waterLRank(spWaterL@data[,"TYPE"],
                                 spWaterL@data[,"NAME"] )      
    spWaterL <- spWaterL[spWaterL$value>=filterVec[4],]
    wLlayer <- shapeToRasterLayer(sxdf=spWaterL,
                                  templateRaster=baseLayer,
                                  maxRasterize=maxRasterize)
    if (!is.finite(wLlayer@data@min)) {
      warning("wLlayer mess-up")
      print(wLlayer@data@min)
    }
  } else {
    wLlayer <- raster::raster(nrows=nrow(baseLayer),
                              ncols=ncol(baseLayer),
                              ext=extent(baseLayer),
                              crs=crs(baseLayer),
                              vals=0)
    print("no water lines to add")
  }
  if (!is.null(spRoads)) {
    print("roads")
    spRoads <- sxdfMask(spRoads,bLrect)
  }
  if (!is.null(spRoads)) {
    print(paste0(nrow(spRoads)," roads to process"))
    spRoads$value <- roadRank(spRoads@data[,"TYPE"],
                              spRoads@data[,"NAME"] )      
    spRoads <- spRoads[spRoads$value>=filterVec[2],]
    rlayer <- shapeToRasterLayer(sxdf=spRoads,
                                 templateRaster=baseLayer,
                                 maxRasterize=maxRasterize)
    if (!is.finite(rlayer@data@min)) {
      warning("rlayer mess-up")
      print(rlayer@data@min)
    }
  } else {
    rlayer <- raster::raster(nrows=nrow(baseLayer),
                             ncols=ncol(baseLayer),
                             ext=extent(baseLayer),
                             crs=crs(baseLayer),
                             vals=0)
    print("no roads to add")
  }
  s <- raster::stack(tlayer,wAlayer,wLlayer,rlayer)
  names(s) <- c("town","waterA","waterL","road")  
  return(s)
}
shapeToRasterLayer <- function(sxdf,templateRaster,
                               maxRasterize=10000,
                               polySimplify=0,polyMethod="vis",
                               polyWeighting=0.85,polySnapInt=0.0001) {
  # return a rasterLayer based on templateRaster, rasterizing sxdf lines/polygons 
  #   using values in sxdf@data[,"rank"]
  zeroRaster <- templateRaster
  raster::values(zeroRaster) <- 0
  retRaster <- zeroRaster
  
  CP <- as(extent(templateRaster), "SpatialPolygons")
  sp::proj4string(CP) <- CRS(sp::proj4string(sxdf))
  sxdf <- sxdfMask(sxdf,CP) 
  
  if (!is.null(sxdf)) {
    sxdf <- sxdf[order(sxdf$value),]  #  sort for velox rasterize
    if (class(sxdf)=="SpatialPolygonsDataFrame") {
      if (polySimplify>0) {
        sxdf <- rmapshaper::ms_simplify(sxdf,keep=polySimplify,
                                        method=polyMethod,
                                        weighting=polyWeighting,
                                        snap=TRUE,snap_interval=polySnapInt) 
        sxdf <- sp::spTransform(sxdf,raster::crs(templateRaster))
      }   
    }
    nloops <- ceiling(nrow(sxdf)/maxRasterize)
    for (i in 1:nloops) {
      if (nloops > 1) print(paste0(i," / ",nloops))
      first <- (i-1)*maxRasterize + 1
      last <- min(i*maxRasterize,nrow(sxdf))
      gc()        #  cleanup, this takes a lot of memory
      rlayer <- velox::velox(zeroRaster)
      print(paste0(round(system.time(
        rlayer$rasterize(spdf=sxdf[first:last,],field="value",background=0)
      )[[3]],digits=2)," rasterizing"))
      if (nloops>1) {
        print(paste0(round(system.time(
          retRaster <- maxLayerValue(retRaster,rlayer$as.RasterLayer(band=1))
        )[[3]],digits=2)," combining"))
      } else {
        retRaster <- rlayer$as.RasterLayer(band=1)
      }
    } 
  }
  return(retRaster)
}
maxLayerValue <- function(rasterLayer1,rasterLayer2) {
  return(max(raster::stack(list(rasterLayer1,rasterLayer2))))
}
townRank <- function(ttype,tname) {
  ttype <- as.vector(ttype)
  # for US TIGER data type is   U pop >= 50k, 2.5k <= C pop < 50k 
  # for Canada type is   B=Metro, K=agglomeration w/tracts, D=agglomeration w/o tracts
  #    (K is dropped already)
  rankT <- rep(1,length(ttype))
  rankT[ttype %in% c("C","D")] <- 3           
  rankT[ttype %in% c("U","B")] <- 5           
  return(as.integer(rankT))
}
roadRank <- function(rtype,rname) {
  rtype <- as.vector(rtype)
  # for US TIGER data type 
  # type is   S1100=Primary        S1200=Secondary      S1400=Local St 
  #           S1500=Vehicular Trl  S1630=Ramp           S1640=Service Drive 
  #           S1710=Walkway        S1720=Stairway       S1730=Alley
  #           S1740=Private Rd     S1750=Internal       S1780=Pkgng Lot 
  #           S1820=Bike Path      S1830=Bridle Trl
  # for Canada type is from(CLASS)
  #   10 Highway, 11 Expressway, 12 Primary highway, 13 Secondary highway
  #   20 Road, 21 Arterial, 22 Collector, 23 Local, 24 Alley/Lane/Utility
  #   25 Connector/Ramp, 26 Reserve/Trail, 27 Rapid transit
  #   80 - bridge/tunnel 
  rankR <- rep(1,length(rtype))                    # anything here gets a 1
  rankR[rtype %in% c("S1630","S1640",
                     "S1730","S1780","S1820",
                     "24","25","26")]         <- 2 # Service Drive, Bike Path
  rankR[rtype %in% c("S1400",
                     "20","21","22","23")]    <- 3 # Local Street
  rankR[rtype %in% c("S1100","10","13","80")] <- 4 # secondary
  rankR[rtype %in% c("S1200","11","12","27")] <- 5 # primary
  return(as.integer(rankR))
}
waterARank <- function(wtype,wname,wsize) {
  wtype <- as.vector(wtype)
  wname <- as.vector(wname)
  wsize <- as.vector(wsize)
  # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
  #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
  # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
  #           CANADA=from canada files
  rankA <- rep(1,length(wtype))                      # anything here gets at least 1
  rankA[wtype %in% c("H2040","H2041","H2060")] <- 2  # res/treatmentpond/pit/quarry
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA")]                <- 3  # lake/pond/swamp/stream
  
  rankA[wtype %in% c("H2040","H2041","H2060",
                     "H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          wsize >= 1000]                        <- 4 # any of the previous that are not small
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          !is.na(wname)]                        <- 5 # lake/pond/swamp/stream named
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          wsize >= 4000]                        <- 6 # lake/pond/swamp/stream big
  rankA[wtype %in% c("H2053","H2051")]          <- 7 # Ocean/Bay/Est/Sound
  rankA[wtype == "H2081"]                       <- 8 # glacier
  return(as.integer(rankA))
}
waterLRank <- function(wtype,wname) {
  wtype <- as.vector(wtype)
  wname <- as.vector(wname)
  #  classify the linear water by ranking and sort it so higher rank overwrites lower
  # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
  #   H1100 is unspecified, gets  set to 1
  rankL <- rep(1,length(wtype))               # anything here gets a 1
  rankL[wtype == "H3020"] <- 2                # canal/ditch
  rankL[wtype == "H3013"] <- 3                # H3013 braided stream
  rankL[wtype %in% c("H3010","CANADA")] <- 4  # H3010 stream/river
  rankL[wtype %in% c("H3010","CANADA") & 
          !is.na(wname)] <- 5                 # H3010 + name not missing
  rankL[wtype %in% c("H3010","CANADA") &
          grepl("RIV",toupper(wname))] <- 6   # H3010 + name contains "RIV"
  return(as.integer(rankL))
}
loadShapeFiles <- function(USStatevec,CAProvincevec,mapshape,
                           shapefileDir,writeShapefiles,
                           shapefileSource="Shapefiles",
                           includeAllRoads=FALSE) {
  workProj4 <- raster::crs(mapshape)
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  if (!is.null(USStatevec)) {
    tmp <- USFeatures(USStatevec=USStatevec,
                      workProj4=workProj4,
                      shapefileSource=shapefileSource,
                      shapefileDir=shapefileDir,
                      writeShapefiles=writeShapefiles,
                      includeAllRoads=includeAllRoads)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  if (!is.null(CAProvincevec)) {
    tmp <- CAFeatures(CAProvincevec,workProj4,shapefileDir)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  #  Don't need to check boundaries - 
  #   saved chunks were appropriately masked/filtered and 
  #   TIGER is masked if needed when fetched
  #   canada data was categorized when shapefiles were created
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
USFeatures <- function(USStatevec,workProj4,
                       shapefileDir,writeShapefiles,
                       shapefileSource="TIGER",
                       includeAllRoads=FALSE) {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  for (st in USStatevec) {
    if (shapefileSource != "Shapefiles" | 
        !file.exists(paste0(shapefileDir,"/",st,"Town.shp"))) {
      tmp <- USTigerFeatures(st,workProj4=workProj4,
                             shapefileDir=shapefileDir,
                             writeShapefiles=writeShapefiles,
                             includeAllRoads=includeAllRoads)
    } else {
      tmp <- readShapeFiles(st,shapefileDir,workProj4)
    }
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }   
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
USTigerFeatures <- function(st,workProj4,
                            shapefileDir,writeShapefiles=TRUE,
                            includeAllRoads=FALSE) {
  stMask <- sp::spTransform(tigris::counties(st),workProj4)
    
  spTown <- sp::spTransform(tigris::urban_areas(),workProj4)
  spTown <- spTown[,(names(spTown) %in% c("NAME10","UATYP10"))]
  colnames(spTown@data) <- c("NAME","TYPE")
  spTown <- sxdfMask(spTown,stMask,keepTouch=TRUE) #keep if town touches the state
  # type is   U pop >= 50k, 2.5k <= C pop < 50k   
    
  # spRoads <- sp::spTransform(tigris::primary_secondary_roads(st),workProj4) # SpatialPolygonsDF
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
                           sp::spTransform(tigris::area_water(st,c),workProj4)) 
      #spatialLines dataframe
      tmpL <- rbind_NULLok(tmpL,
                           sp::spTransform(tigris::linear_water(st,c),workProj4)) 
      tmpR <- rbind_NULLok(tmpR,
                           sp::spTransform(tigris::roads(st,c),workProj4)) 
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
CAFeatures <- function(CAProvincevec,workProj4,shapefileDir) {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  ##  files are already filtered for city/road/polygonwater importance
  for (pr in CAProvincevec) {
    tmp <- readShapeFiles(pr,shapefileDir,workProj4)
    spTown <- rbind_NULLok(spTown,tmp[["spTown"]])
    spRoads <- rbind_NULLok(spRoads,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spWaterA,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spWaterL,tmp[["spWaterL"]])
  }
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
