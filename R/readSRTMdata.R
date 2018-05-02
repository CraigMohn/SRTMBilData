loadMapElevData <- function(mapcrop,mapDataDir,resstr) {
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
      } else if (rgeos::gIntersects(mapcrop, ei)) {
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
