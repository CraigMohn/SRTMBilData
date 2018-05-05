buildFeatureStack <- function(baseLayer,mapshape,
                              spTown,spWaterA,spWaterL,spRoads,
                              filterVec=c(1,1,1,1),  #towns,roads,waterA,waterL
                              polySimplify=0,polyMethod="vis",
                              polyWeighting=0.85,polySnapInt=0.0001  ) {
  #  baseLayer a rasterLayer
  #  returns a rasterStack  which has 4 layers which can be stored as ints
  bLrect <- as(raster::extent(baseLayer), "SpatialPolygons")
  sp::proj4string(bLrect) <- sp::proj4string(baseLayer)
  
  if (!is.null(spTown)) {
    print("towns")
    tspTown <- sxdfMask(spTown, bLrect)
    if (!is.null(tspTown)) {
      print(paste0(nrow(tspTown)," towns to process"))
      tspTown$value <- townRank(tspTown@data[,"TYPE"],
                               tspTown@data[,"NAME"])
      tspTown <- tspTown[tspTown$value>=filterVec[1],]
      tlayer <- shapeToRasterLayer(sxdf=tspTown,
                                   templateRaster=baseLayer,
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
  }
  if (!is.null(spWaterA)) {
    print("water polygons")
    tspWaterA <- sxdfMask(spWaterA,bLrect)
    if (!is.null(tspWaterA)) {
      print(paste0(nrow(tspWaterA)," water areas to process"))
      tspWaterA$value <- waterARank(tspWaterA@data[,"TYPE"],
                                    tspWaterA@data[,"NAME"],
                                    tspWaterA@data[,"size"] )      
      tspWaterA <- tspWaterA[tspWaterA$value>=filterVec[3],]
      wAlayer <- shapeToRasterLayer(sxdf=tspWaterA,
                                    templateRaster=baseLayer,
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
  }
  if (!is.null(spWaterL)) {
    print("water lines")
    tspWaterL <- sxdfMask(spWaterL,bLrect)
    if (!is.null(tspWaterL)) {
      print(paste0(nrow(tspWaterL)," water lines to process"))
      tspWaterL$value <- waterLRank(tspWaterL@data[,"TYPE"],
                                    tspWaterL@data[,"NAME"] )      
      tspWaterL <- tspWaterL[tspWaterL$value>=filterVec[4],]
      wLlayer <- shapeToRasterLayer(sxdf=tspWaterL,
                                    templateRaster=baseLayer)
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
  }
  if (!is.null(spRoads)) {
    print("roads")
    tspRoads <- sxdfMask(spRoads,bLrect)
    if (!is.null(tspRoads)) {
      print(paste0(nrow(tspRoads)," roads to process"))
      tspRoads$value <- roadRank(tspRoads@data[,"TYPE"],
                                tspRoads@data[,"NAME"] )      
      tspRoads <- tspRoads[tspRoads$value>=filterVec[2],]
      rlayer <- shapeToRasterLayer(sxdf=tspRoads,
                                   templateRaster=baseLayer)
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
  }
  s <- raster::stack(tlayer,wAlayer,wLlayer,rlayer)
  names(s) <- c("town","waterA","waterL","road")  
  return(s)
}
shapeToRasterLayer <- function(sxdf,templateRaster,
                               polySimplify=0,polyMethod="vis",
                               polyWeighting=0.85,polySnapInt=0.0001,
                               maxShapes=200000) {
  # return a rasterLayer based on templateRaster, rasterizing sxdf lines/polygons 
  #   using values in sxdf@data[,"rank"]
  zeroRaster <- templateRaster
  raster::values(zeroRaster) <- 0
  rlayer <- velox::velox(zeroRaster)
  #sxdf <- raster::crop(sxdf,raster::extent(templateRaster))
  #sxdf <- raster::intersect(sxdf,
                        #    as(raster::extent(templateRaster),'SpatialPolygons'))
  CP <- as(extent(templateRaster), "SpatialPolygons")
  sp::proj4string(CP) <- CRS(sp::proj4string(templateRaster))
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
    nchunks <- ceiling(nrow(sxdf)/maxShapes)
    for ( i in 1:nchunks ) {
      gc()        #  cleanup, this takes a lot of memory
      upperlimit <- min(i*maxShapes,nrow(sxdf))
      tsxdf <- sxdf[((i-1)*maxShapes+1):upperlimit,]
      print(system.time(
        rlayer$rasterize(spdf=tsxdf,field="value",background=0)
      )[[3]])
    }
  } 
  return(rlayer$as.RasterLayer(band=1))
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
  # for US TIGER data type is S1100=secondary   S1200=Primary
  # for Canada type is from(RANK) 1=transcanada, 2=national, 3=major, 4=secondary hwy
  rankR <- rep(1,length(rtype))                       # anything here gets a 1
  rankR[rtype %in% c("S1400","S1640","S1820")] <- 2   # Local Street,Service Drive, Bike Path
  rankR[rtype %in% c("S1100","3","4")] <- 3           # secondary
  rankR[rtype %in% c("S1200","1","2")] <- 5             # primary
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
loadShapeFiles <- function(USStatevec,CAProvincevec,mapcrop,
                           shapefileDir,writeShapefiles,
                           shapefileSource="Shapefiles",
                           includeAllRoads=FALSE) {
  workProj4 <- raster::crs(mapcrop)
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
    spRoads <- rbind_NULLok(spTown,tmp[["spRoads"]])
    spWaterA <- rbind_NULLok(spTown,tmp[["spWaterA"]])
    spWaterL <- rbind_NULLok(spTown,tmp[["spWaterL"]])
  }
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}

