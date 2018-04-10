yRatio <- function(rrr) {
  xmin <- rrr@extent@xmin
  xmax <- rrr@extent@xmax
  ymin <- rrr@extent@ymin
  ymax <- rrr@extent@ymax
  return(yRatioPts(xmin,xmax,ymin,ymax))
}
yRatioPts <- function(xmin,xmax,ymin,ymax) {
  width <-
    (raster::pointDistance(cbind(xmin,ymin),cbind(xmax,ymin),lonlat=TRUE) +
       raster::pointDistance(cbind(xmin,ymax),cbind(xmax,ymax),lonlat=TRUE)) / 2
  height <-
    (raster::pointDistance(cbind(xmin,ymin),cbind(xmin,ymax),lonlat=TRUE) +
       raster::pointDistance(cbind(xmax,ymin),cbind(xmax,ymax),lonlat=TRUE)) / 2
  return(height/width)
}
pan3d <- function(button) {
  start <- list()
  begin <- function(x, y) {
    start$userMatrix <<- rgl::par3d("userMatrix")
    start$viewport <<- rgl::par3d("viewport")
    start$scale <<- rgl::par3d("scale")
    start$projection <<- rgl::rgl.projection()
    start$pos <<- rgl::rgl.window2user( x/start$viewport[3], 
                                        1 - y/start$viewport[4], 
                                        0.5,
                                        projection = start$projection)
  }
  update <- function(x, y) {
    xlat <- (rgl::rgl.window2user( x/start$viewport[3], 
                                   1 - y/start$viewport[4], 
                                   0.5,
                                   projection = start$projection) - start$pos)*start$scale
    mouseMatrix <- rgl::translationMatrix(xlat[1], xlat[2], xlat[3])
    rgl::par3d(userMatrix = start$userMatrix %*% t(mouseMatrix) )
  }
  rgl::rgl.setMouseCallbacks(button, begin, update)
  cat("Callbacks set on button", button, "of rgl device", rgl.cur(), "\n")
}

addmapfiles <- function(filenames,lonChar,minLon,maxLon,latChar,minLat,maxLat,
                        resstr="_1arc_v3_bil") {
  #  first pass assume N and W quartersphere
  if (minLat <= maxLat) {
    for (lat in seq(minLat,maxLat)) {
      if (minLon <= maxLon) {
        for (lon in seq(minLon,maxLon))
          filenames <- c(filenames,
                         paste0(latChar,sprintf("%02d",lat),"_",
                                lonChar,sprintf("%03d",lon),
                                resstr,".zip"))
      }
    }
  }
  return(filenames)  
}
bufferUnion <- function(spObj,mapbuffer,mapunion,
                        outCRS="+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                        bufferCRS="+init=epsg:3857",
                        capStyle="FLAT",
                        joinStyle="BEVEL",
                        simplifytol=0) {
  if (mapbuffer > 0) {
    spObj <- rgeos::gBuffer(sp::spTransform( spObj, CRS( bufferCRS ) ),
                            width=1.5*mapbuffer,
                            capStyle=capStyle)
    if (simplifytol > 0) spObj <- rgeos::gSimplify(spObj,tol=simplifytol)
    spObj <- rgeos::gBuffer(spObj,width=-mapbuffer/2)
  }    
  spObj <- sp::spTransform( spObj, CRS( outCRS ) ) 
  if (is.null(mapunion)) {
    return(spObj)
  } else {
    return(rgeos::gUnaryUnion(raster::union(mapunion,spObj)))
  } 
}

rbind_NULLok <- function(a,b) {
  if (is.null(a)) {
    return(b)
  } else {
    return(rbind(a,b))
  }
}

doubleExtent <- function(objectWithExtent) {
  bb <- sp::bbox(objectWithExtent)
  return(extent(max(-180,(3*bb[1,1]-bb[1,2])/2),
                min(180,(-bb[1,1]+3*bb[1,2])/2),
                max(-90,(3*bb[2,1]-bb[2,2])/2),
                min(90,(-bb[2,1]+3*bb[2,2])/2)))
}

loadStateElevData <- function(USStatevec,CAProvincevec) {
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  spTown <- NULL
  j <- 1
  r.list <- list()
  for (st in c(USStatevec,CAProvincevec)) {
    fvec <- list.files(path=paste0(datadir,"/rasterfiles"),
                       pattern=paste0(st,"elevs[0-9]{,2}.grd"))
    for (fn in fvec) {
      print(paste0("loading ",fn))
      elevations <- raster(paste0(datadir,"/rasterfiles/",fn))
      print(elevations)
      r.list[[j]] <- elevations
      j <- j + 1
    }
    sdfstr <- paste0(datadir,"/shapefiles/",st)
    tRoads <- sp::spTransform(raster::shapefile(paste0(sdfstr,"Roads.shp")),
                              workProj4)
    spRoads <- rbind_NULLok(spRoads,tRoads)
    tWaterA <- sp::spTransform(raster::shapefile(paste0(sdfstr,"WaterA.shp")),
                               workProj4)
    spWaterA <- rbind_NULLok(spWaterA,tWaterA)
    tWaterL <- sp::spTransform(raster::shapefile(paste0(sdfstr,"WaterL.shp")),
                               workProj4)
    spWaterL <- rbind_NULLok(spWaterL,tWaterL)
    tTown <- sp::spTransform(raster::shapefile(paste0(sdfstr,"Town.shp")),
                             workProj4)
    spTown <- rbind_NULLok(spTown,tTown)
  }
  print("calling merge")
  if (j > 2) {
    m.sub <- do.call(merge, r.list)
  } else {
    m.sub <- elevations
  }
  #m.sub <- raster::mask(raster::crop(m.sub, extent(mapcrop),snap="near"),
  #                      mapcrop)
  return(list(elevraster=m.sub,spTown=spTown,spRoads=spRoads,
              spWaterA=spWaterA,spWaterL=spWaterL))
}

loadMapElevData <- function(mapcrop) {
  m.sub <- NULL
  j <- 1
  r.list <- list()
  mapextent <- raster::extent(mapcrop)
  print(mapextent)
  ELonMin <- floor(max(mapextent@xmin,0))
  ELonMax <- floor(max(mapextent@xmax,0))
  WLonMin <- ceiling(max(-mapextent@xmax,0))
  WLonMax <- ceiling(max(-mapextent@xmin,0))
  NLatMin <- floor(max(mapextent@ymin,0))
  NLatMax <- floor(max(mapextent@ymax,0))
  SLatMin <- ceiling(max(-mapextent@ymax,0))
  SLatMax <- ceiling(max(-mapextent@ymin,0))
  
  fn <- NULL
  if (NLatMax > 0) {
    if (WLonMax > 0) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"n",NLatMin,NLatMax,resstr)
    }
    if (ELonMax > 0) {
      fn <- addmapfiles(fn,"e",ELonMin,ELonMax,"n",NLatMin,NLatMax,resstr)
    }
  }
  if (SLatMax > 0) {
    if (WLonMax > 0) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"s",SLatMin,SLatMax,resstr)
    }
    if (ELonMax > 0) {
      fn <- addmapfiles(fn,"e",ELonMin,ELonMax,"s",SLatMin,SLatMax,resstr)
    }
  }
  tempd <- tempdir()
  rn <- gsub("_bil","",tools::file_path_sans_ext(fn))
  unfoundfn <- NULL
  firstRes <- NULL
  firstProj <- NULL
  firstOrigin <- NULL
  for (i in 1:length(fn)) {
    if (file.exists(paste0(mapDataDir,"/",fn[[i]]))) {
      cat("\n")
      print(paste0(mapDataDir,"/",fn[[i]]))
      unzip(paste0(mapDataDir,"/",fn[[i]]),exdir=tempd)
      tmp <- raster(paste0(tempd,"/",rn[[i]],".bil"))
      if (is.null(firstOrigin)) {
        print(tmp)
        firstOrigin <- raster::origin(tmp)
        print(firstOrigin)
        firstXWide <- tmp@extent@xmax - tmp@extent@xmin
        firstYWide <- tmp@extent@ymax - tmp@extent@ymin
        firstRows <- nrow(tmp)
        firstCols <- ncol(tmp)
        firstProj <- raster::projection(tmp)
        firstRes <- raster::res(tmp)
      } else {
        if (raster::projection(tmp)!=firstProj) 
          warning("projection mismatch - ",raster::projection(tmp))
      }
      if (!identical(firstRes,raster::res(tmp))) {
        print("resolution differs from the first tile - resampling")
        print(raster::res(tmp))
        ## what raster do we want? - clone first in dims, extent size and offset 
        xwide <- tmp@extent@xmax - tmp@extent@xmin
        ywide <- tmp@extent@ymax - tmp@extent@ymin
        llx <- tmp@extent@xmin - (firstXWide-xwide)/2
        lly <- tmp@extent@ymin - (firstYWide-ywide)/2
        newraster <- tmp
        raster::ncol(newraster) <- firstCols
        raster::nrow(newraster) <- firstRows
        newraster <- raster::setExtent(newraster,
                                       extent(llx,llx+firstXWide,
                                              lly,lly+firstYWide))
        raster::origin(newraster) <- firstOrigin
        print(tmp)  
        tmp <- raster::resample(tmp, newraster)
        print(tmp)
      } else {
        if (max(abs(firstOrigin-raster::origin(tmp))) > 0.0000001)  
          warning("origin mismatch - ",raster::origin(tmp)," ",firstOrigin)
      }
      xmin <- tmp@extent@xmin
      xmax <- tmp@extent@xmax
      ymin <- tmp@extent@ymin
      ymax <- tmp@extent@ymax
      pgon <- sp::Polygon(cbind(c(xmin,xmax,xmax,xmin,xmin),
                                c(ymin,ymin,ymax,ymax,ymin)))
      ei <- sp::SpatialPolygons(list(Polygons(list(pgon), ID = "1deg")),
                                proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
      if (rgeos::gContainsProperly(mapcrop, ei)) {
        print("interior - not masked")
      } else if (gIntersects(mapcrop, ei)) {
        print(paste0("boundary - masking time = ",system.time(
          tmp <- raster::mask(raster::crop(tmp, extent(mapcrop),snap="near"),
                              mapcrop)
        ))[[3]])
        print(origin(tmp))
      } else {
        print("exterior - not used")
        tmp <- NULL
      }
      if (!is.null(tmp)) {
        r.list[[j]] <- tmp
        j <- j + 1
      }
    } else {
      print(paste0(mapDataDir,"/",fn[[i]]," does not exist, ignored"))
      unfoundfn <- c(unfoundfn,fn[[i]])
    }
  }
  print(warnings())
  print("calling merge")
  if (j > 2) {
    #m.sub <- do.call(merge, r.list))
    m.sub <- do.call(merge, c(r.list,list(tolerance=0.1)))
  } else {
    m.sub <- r.list[[1]]
  }
  return(m.sub)
}



USFeatures <- function(USStatevec,workProj4,
                       writeShapefiles=FALSE,shapefiledir) {
  spTown <- NULL
  spRoads <- NULL
  spWaterA <- NULL
  spWaterL <- NULL
  
  spTownUS <- sp::spTransform(tigris::urban_areas(),workProj4) # SpatialPolygonsDF
  spTownUS <- spTownUS[,(names(spTownUS) %in% c("NAME10","UATYP10"))]
  colnames(spTownUS@data) <- c("NAME","TYPE") # type is   A pop >= 50k, 2.5k <= U pop < 50k
  for (st in USStatevec) {
    stMask <- sp::spTransform(tigris::counties(st),sp::CRS(workProj4))
    tmpT <- sxdfMask(spTownUS,stMask,keepTouch=TRUE) #keep if town touches the state
    # type is   A pop >= 50k, 2.5k <= U pop < 50k    
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
    tmpA <- tmpA[,(names(tmpA) %in% c("FULLNAME","MTFCC"))]
    colnames(tmpA@data) <- c("NAME","TYPE")  
    # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
    #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
    tmpL <- tmpL[,(names(tmpL) %in% c("FULLNAME","MTFCC"))]
    colnames(tmpL@data) <- c("NAME","TYPE")  
    # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
    tmpL <- tmpL[tmpL@data[,"TYPE"]=="H3010",]  # keep only rivers and streams H3010

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
sxdfMask <- function(sxdf,poly,keepTouch=FALSE) {
  if (is.null(sxdf) | is.null(poly)) return(NULL) 
  
  tmpdata <- sxdf@data
  
  tmp.1 <- rgeos::gIntersects(sxdf, poly, byid=TRUE)
  tmp.2 <- as.logical(apply(tmp.1, 2, function(x) {sum(x)} ))
  
  #  keep only intersecting sp objects in dataframe
  tmpgeo <- sxdf[tmp.2,]
  tmpdata <- tmpdata[tmp.2,]
  if (!keepTouch) tmpgeo <- rgeos::gIntersection(tmpgeo, poly, byid=TRUE)
  row.names(tmpgeo) <- row.names(tmpdata)
  if (class(sxdf)=="SpatialLinesDataFrame") {
    return(sp::SpatialLinesDataFrame(tmpgeo,data=tmpdata))
  } else if (class(sxdf)=="SpatialPolygonsDataFrame") {
    return(sp::SpatialPolygonsDataFrame(tmpgeo,data=tmpdata))
  } else {
    stop(" sxdfMask expects a spatialLinesDataFrame or a spatialPolygonsDataFrame")
  }
}

shapes_for_states <- function(statevec,
                              workProj4="+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs",
                              shapefiledir="c:/bda/shapefiles") {
  tmp <- USFeatures(statevec,workProj4,
                         writeShapefiles=TRUE,shapefiledir)
  return("done")
}

