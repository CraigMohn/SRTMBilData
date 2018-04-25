USFeatures <- function(USStatevec,workProj4,
                       writeShapefiles=FALSE,shapefiledir) {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  
  spTownUS <- sp::spTransform(tigris::urban_areas(),workProj4) # SpatialPolygonsDF
  spTownUS <- spTownUS[,(names(spTownUS) %in% c("NAME10","UATYP10"))]
  colnames(spTownUS@data) <- c("NAME","TYPE")
  for (st in USStatevec) {
    stMask <- sp::spTransform(tigris::counties(st),sp::CRS(workProj4))
    tmpT <- sxdfMask(spTownUS,stMask,keepTouch=TRUE) #keep if town touches the state
    # type is   U pop >= 50k, 2.5k <= C pop < 50k    
    spTown <- rbind_NULLok(spTown,tmpT)
    plot(tmpT)    
    tmpR <- sp::spTransform(tigris::primary_secondary_roads(st),workProj4) # SpatialPolygonsDF
    tmpR <- tmpR[,(names(tmpR) %in% c("FULLNAME","MTFCC"))]
    colnames(tmpR@data) <- c("NAME","TYPE")  
    # type is   S1100=secondary   S1200=Primary
    spRoads <- rbind_NULLok(spRoads,tmpR)
    plot(tmpR)    
    
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
      }
    }  
    tmpA$size <- as.numeric(tmpA$AWATER) + as.numeric(tmpA$ALAND)
    tmpA <- tmpA[,(names(tmpA) %in% c("FULLNAME","MTFCC","size"))]
    colnames(tmpA@data) <- c("NAME","TYPE","size")  
    # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
    #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
    tmpL <- tmpL[,(names(tmpL) %in% c("FULLNAME","MTFCC"))]
    colnames(tmpL@data) <- c("NAME","TYPE")  
    # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
    #tmpL <- tmpL[tmpL@data[,"TYPE"]=="H3010",]  # keep only rivers and streams H3010

    if (writeShapefiles) writeFeatures(st,shapefiledir=shapefiledir,
                                       spTown=tmpT,spRoads=tmpR,
                                       spWaterA=tmpA,spWaterL=tmpL)
plot(tmpA)
plot(tmpL)
    spWaterA <- rbind_NULLok(spWaterA,tmpA)
    spWaterL <- rbind_NULLok(spWaterL,tmpL)
  }
  
  if (!is.null(spTown)) plot(spTown)
  if (!is.null(spRoads)) plot(spRoads)
  if (!is.null(spWaterA)) plot(spWaterA)
  if (!is.null(spWaterL)) plot(spWaterL)

  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
CAFeatures <- function(CAProvincevec,workProj4) {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  ##  files are already filtered for city/road/polygonwater importance
  for (pr in CAProvincevec) {
    print(paste0("loading features for ",pr))
    if (showCities) {
      tmp <- raster::shapefile(paste0(datadir,"/shapefiles/",pr,
                                      "Town.shp")) # SpatialPolygonsDF
      tmp <- sp::spTransform(tmp,workProj4)
      spTown <- rbind_NULLok(spTown,tmp)
     }
    if (showRoads) {
      tmp <- raster::shapefile(paste0(datadir,"/shapefiles/",pr,
                                      "Roads.shp")) # SpatialLinesDF
      tmp <- sp::spTransform(tmp,workProj4)
      spRoads <- rbind_NULLok(spRoads,tmp)
    }
    if (showWater) {
      tmp <- raster::shapefile(paste0(datadir,"/shapefiles/",pr,
                                      "WaterA")) # SpatialPolygonsDF
      tmp <- sp::spTransform(tmp,workProj4)
      spWaterA <- rbind_NULLok(spWaterA,tmp)
      tmp <- raster::shapefile(paste0(datadir,"/shapefiles/",pr,
                                      "WaterL")) # SpatialLinesDF
      tmp <- sp::spTransform(tmp,workProj4)
      spWaterL <- rbind_NULLok(spWaterL,tmp)
    }
  }
  if (!is.null(spTown)) plot(spTown, col="PeachPuff")
  if (!is.null(spRoads)) plot(spRoads)
  if (!is.null(spWaterA)) plot(spWaterA)
  if (!is.null(spWaterL)) plot(spWaterL)
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}
writeFeatures <- function(stname,shapefiledir,spTown,spRoads,spWaterA,spWaterL) {
  print(paste0("writing feature shapefiles for ",stname))
  raster::shapefile(spRoads,filename=paste0(shapefiledir,"/",stname,"Roads.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterA,filename=paste0(shapefiledir,"/",stname,"WaterA.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterL,filename=paste0(shapefiledir,"/",stname,"WaterL.shp"),
                    overwrite=TRUE)
  raster::shapefile(spTown,filename=paste0(shapefiledir,"/",stname,"Town.shp"),
                    overwrite=TRUE)
  return(NULL)
}
filterWaterA <- function(waterLineDF,level="named") {
  if (level == "named") {
    tmpkeep <- !is.na(waterLineDF@data[,"NAME"])
    return(waterLineDF[tmpkeep,])
  } else if (level=="lakes") {
    tmpkeep <- waterLineDF@data[,"TYPE"] %in% c("H2030")
    return(waterLineDF[tmpkeep,])
  } else if (level=="noarea") {
    return(NULL)
  } else {
    return(waterLineDF)
  } 
}
filterWaterL <- function(waterLineDF,level="area") {
  waterLineDF <- waterLineDF[waterLineDF$TYPE %in% c("H3010",""),]
  if (is.null(waterLineDF)) return(NULL)
  if (level=="area" | level == "lakes") {
    return(NULL)
  } else if (level=="all") {
    return(waterLineDF)
  } else {
    tmpname <- waterLineDF@data[,"NAME"]
    if ((level == "rivers") | (level == "noarea")) {
      tmpkeep <- grepl("Riv",tmpname)
    } else if (level == "riversplus") {
      tmpkeep <- !is.na(tmpname) & 
        !grepl("Ditch",tmpname) &
        !grepl("Gulch",tmpname) &
        !grepl("Cyn",tmpname) &
        !grepl("Crk",tmpname) &
        !grepl("Tributary",tmpname) &
        !grepl("Watercourse",tmpname) &
        !grepl("Strm",tmpname)
    } else if (level == "named") {
      tmpkeep <- !is.na(tmpname)
    } else {
      stop("bad value for level in filterWaterL")
    }
    return(waterLineDF[tmpkeep,])
  }
}

buildFeatureStack <- function(baseLayer,mapshape,
                              spTown,spWaterA,spWaterL,spRoads,
                              fastAreas=TRUE,
                              polySimplify=0,polyMethod="vis",
                              polyWeighting=0.85,polySnapInt=0.0001  ) {
  #  baseRaster a raster, rasterStack or rasterBrick
  #  returns a rasterStack  which has 4 layers which can be stored as ints
  #  - to do - codes for canada
  
  tlayer <- raster::raster(nrows=nrow(baseLayer),
                           ncols=ncol(baseLayer),
                           ext=extent(baseLayer),
                           crs=crs(baseLayer),
                           vals=0)
  if (!is.null(spTown)) {
    print("towns")
    tspTown <- raster::intersect(spTown, tlayer)
    if (!is.null(tspTown)) {
      # type is   U pop >= 50k, 2.5k <= C pop < 50k 
      townlevel = ifelse(tspTown@data[,"TYPE"] %in% c("U","B"),5,3)
      print(paste0(nrow(tspTown)," towns to process"))
      if (fastAreas) {
        temp <- rgeos::gUnaryUnion(tspTown)
        tlayer <- raster::rasterize(temp,tlayer,
                                    field=townlevel,update=TRUE)
      } else {
        print(paste0(nrow(tspTown)," towns to process"))
        for (i in 1:nrow(tspTown)) {
          cat("\r",i,"  ",tspTown@data[i,"NAME"],
                     "                                    ")
          tlayer <- raster::rasterize(tspTown[i,],tlayer,
                                    field=townlevel,update=TRUE)
        }
        cat("\n")
      }
      plot(tspTown,col=citycolor)
      if (!is.finite(tlayer@data@min)) {
        warning("tlayer mess-up")
        print(tlayer@data@min)
      }
      print(paste0(nrow(tspTown)," towns rasterized"))
    } else {
      print("no towns to add")
    }
  }
  
  wAlayer <- raster::raster(nrows=nrow(baseLayer),
                           ncols=ncol(baseLayer),
                           ext=extent(baseLayer),
                           crs=crs(baseLayer),
                           vals=0)
  if (!is.null(spWaterA)) {
    print("water polygons")
    tspWaterA <- raster::intersect(spWaterA, wAlayer)
    if (!is.null(tspWaterA)) {
      #  classify the area water by ranking and sort it so higher rank overwrites lower
      ttype <- tspWaterA@data[,"TYPE"]
      # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
      #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
      #  need field for area of water body
      rankA <- rep(1,length(ttype))                      # anything here gets at least 1
      rankA[ttype == "H2081"] <- 6                       # glacier
      rankA[ttype %in% c("H2053","H2051")] <- 5          # Ocean/Bay/Est/Sound
      rankA[ttype %in% c("H2025","H2030")] <- 4          # lake/pond/swamp
      rankA[ttype %in% c("H2025","H2030") &
            tspWaterA@data[,"size"] < 4000] <- 3         # lake/pond/swamp(<4000 meters2)
      rankA[ttype %in% c("H2040","H2041","H2060")] <- 2  # res/treatmentpond/pit/quarry
      tspWaterA$rank <- rankA      
      tspWaterA <- tspWaterA[order(tspWaterA$rank),]     ### sort it
      print(paste0(nrow(tspWaterA)," water areas to process"))
      if (fastAreas) {
        for (i in  sort(unique(rankA))) {
          print(paste0("rank=",i,"  ", 
                       nrow(tspWaterA[tspWaterA@data[,"rank"]==i,]),
                       " areas to process"))
          temp <- rgeos::gUnaryUnion(tspWaterA[tspWaterA@data[,"rank"]==i,])
          wAlayer <- raster::rasterize(temp,wAlayer,
                                       field=i,update=TRUE)
        }
        nWaterA <- nrow(tspWaterA)
      } else {
        nWaterA <- 0
        for (i in 1:nrow(tspWaterA)) {
          if (polySimplify == 0) {
            temp <- rgeos::gSimplify(tspWaterA[i,],tol=0,topologyPreserve=TRUE)
          } else {
            temp <- rmapshaper::ms_simplify(tspWaterA[i,],keep=polySimplify,
                                            method=polyMethod,
                                            weighting=polyWeighting,
                                            snap=TRUE,snap_interval=polySnapInt) 
            temp <- sp::spTransform(temp,raster::crs(wAlayer))
          }
          if (rgeos::gIntersects(temp,mapshape)) {
            cat("\r",i,"  ",tspWaterA@data[i,"NAME"],
                       "                                 ")
            temp <- rgeos::gIntersection(temp,mapshape,drop_lower_td=TRUE)
            if (!is.null(temp)) {
              wAlayer <- raster::rasterize(temp,wAlayer,
                                           field=tspWaterA@data[i,"rank"],update=TRUE)
              nWaterA <- nWaterA + 1
            }
          }
        }
        cat("\n")
      }
      plot(tspWaterA, col="blue")
      if (!is.finite(wAlayer@data@min)) {
        warning("wAlayer mess-up")
        print(wAlayer@data@min)
      }
      print(paste0(nWaterA," water polygons drawn"))
    } else {
      print("no water polygons to add")
    }
  }
  
  wLlayer <- raster::raster(nrows=nrow(baseLayer),
                            ncols=ncol(baseLayer),
                            ext=extent(baseLayer),
                            crs=crs(baseLayer),
                            vals=0)
  if (!is.null(spWaterL)) {
    print(paste0("water lines - ",nrow(spWaterL)))
    #  classify the linear water by ranking and sort it so higher rank overwrites lower
    ttype <- spWaterL@data[,"TYPE"]
    tname <- spWaterL@data[,"NAME"]
    # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
    rankL <- rep(1,length(ttype))                      # anything here gets a 1
    rankL[ttype == "H3020"] <- 2                       # canal/ditch
    rankL[ttype == "H3013"] <- 3                       # H3013 braided stream
    rankL[ttype == "H3010"] <- 4                       # H3010 stream/river
    rankL[ttype == "H3010" &
            grepl("RIV",toupper(tname))] <- 6            # H3010 + name contains "RIV"
    spWaterL$rank <- rankL      

    for (i in  sort(unique(rankL))) {
      tspWaterL <- spWaterL[spWaterL$rank == i,]
      nWaterLi <- nrow(tspWaterL)
      nWaterLii <- 0
      for (j in 1:ceiling(nWaterLi/20000)) { 
        keepidx <- seq((j-1)*20000+1,min(nWaterLi,j*20000),1)
        ttspWaterL <- tspWaterL[keepidx,]
        ttspWaterL <- raster::crop(ttspWaterL,wLlayer)
        if (!is.null(ttspWaterL)) {
          nWaterLii <- nWaterLii + nrow(ttspWaterL)
          #temp <- rgeos::gLineMerge(ttspWaterL)
          temp <- ttspWaterL
print(system.time(
          wLlayer <- raster::rasterize(temp,wLlayer,
                                     field=i,update=TRUE)
)[3])
        }
      }
      print(paste0("rank ",i,"  added ",nWaterLii,
                   " linear water features out of list of ", nWaterLi))
    }
    plot(tspWaterL, col="blue")
    if (!is.finite(wLlayer@data@min)) {
      warning("wLlayer mess-up")
      print(wLlayer@data@min)
    }
  } else {
    print("no water lines to add")
  }
  
  rlayer <- raster::raster(nrows=nrow(baseLayer),
                           ncols=ncol(baseLayer),
                           ext=extent(baseLayer),
                           crs=crs(baseLayer),
                           vals=0)
  if (!is.null(spRoads)) {
    print("roads")
    tspRoads <- raster::crop(spRoads, rlayer) 
    if (!is.null(tspRoads)) {
      nRoads <- nrow(tspRoads)
      # type is   S1100=secondary   S1200=Primary
      ttype <- tspRoads@data[,"TYPE"]
      rankR <- rep(1,length(ttype))                      # anything here gets a 1
      rankR[ttype == "S1100"] <- 3                       # secondary
      rankR[ttype == "S1200"] <- 5                       # primary
      tspRoads$rank <- rankR      
      for (i in sort(unique(rankR))) {
        #temp <- rgeos::gLineMerge(tspRoads[tspRoads@data[,"rank"]==i,])
        temp <- tspRoads[tspRoads@data[,"rank"]==i,]
        rlayer <- raster::rasterize(temp,rlayer,
                                    field=i,update=TRUE)
      }
      plot(tspRoads)
      if (!is.finite(rlayer@data@min)) {
        warning("rlayer mess-up")
        print(rlayer@data@min)
      }
      print(paste0(nRoads," roads drawn"))
    } else {
      print("no roads to add")
    }
  }
  s <- raster::stack(tlayer,wAlayer,wLlayer,rlayer)
  names(s) <- c("town","waterA","waterL","road")  
  return(s)
}
