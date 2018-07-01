loadMapElevData <- function(mapshape,mapDataDir,resstr,
                            noisy=TRUE,silent=FALSE,
                            latLimit=59) {
  m.sub <- NULL
  j <- 1
  r.list <- list()
  mapextent <- raster::extent(mapshape)
  if (!silent) print(mapextent)

  NLatMin <- floor(max(mapextent@ymin,0))
  NLatMax <- floor(max(mapextent@ymax,0))
  SLatMin <- ceiling(max(-mapextent@ymax,0))
  SLatMax <- ceiling(max(-mapextent@ymin,0))

  if (!is.null(mapDataDir))tempd <- tempdir()
  unfoundfn <- NULL
  firstRes <- NULL
  firstProj <- NULL
  firstOrigin <- NULL
  #  loop over N hemisphere from equator, then S from equator (highest res data there)
  latseq <- NULL
  if (NLatMax>0) latseq <- seq(from=NLatMin,to=min(NLatMax,latLimit))
  if (SLatMax>0) latseq <- c(latseq,-seq(from=(SLatMin+1),to=min(SLatMax,latLimit)))
  lonseq <- seq(from=floor(mapextent@xmin),to=floor(mapextent@xmax))

  for (lat in latseq) {
    for (lon in lonseq) {
      pgon <- sp::Polygon(cbind(c(lon,(lon+1),(lon+1),lon,lon),
                                c(lat,lat,(lat+1),(lat+1),lat)))
      ei <- sp::SpatialPolygons(list(Polygons(list(pgon), ID = "1deg")),
                                proj4string=CRS(sp::proj4string(mapshape)))
                                #proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
      if (rgeos::gIntersects(mapshape, ei)) {
        if (is.null(mapDataDir)) {
          #  need to fix to handle no tile returned
          if (!silent) print(paste0(" ",lat," , ",lon))
          tmp <- getData("SRTM",lon=floor(lon)+0.5,lat=floor(lat)+0.5)
        } else {
          fname <- mapfilename(lat,lon,resstr)
          if (file.exists(paste0(mapDataDir,"/",fname))) {
            # cat("\n")
            if (!silent) print(paste0(mapDataDir,"/",fname))
            unzip(paste0(mapDataDir,"/",fname),exdir=tempd)
            rname <- gsub("_bil","",tools::file_path_sans_ext(fname))
            tmp <- raster(paste0(tempd,"/",rname,".bil"))
          } else {
            if (!silent) print(paste0(mapDataDir,"/",fname," does not exist, ignored"))
            tmp <- NULL
          }
        }
        if (!is.null(tmp)) {
          if (is.null(firstOrigin)) {
            firstRaster <- tmp
            if (noisy) print(tmp)
            firstOrigin <- raster::origin(tmp)
            if (noisy) print(firstOrigin)
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
          if (!raster::compareRaster(firstRaster,tmp,
                                     res=TRUE,rowcol=FALSE,
                                     extent=FALSE,orig=FALSE,
                                     stopiffalse=FALSE)) {
            if (noisy) print("resolution differs from the first tile - resampling")
            if (noisy) print(raster::res(tmp))
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
            if (noisy) print(tmp)  
            tmp <- raster::resample(tmp, newraster)
            if (noisy) print(tmp)
          } else {
            if (max(abs(firstOrigin-raster::origin(tmp))) > 0.0000001)  
              warning("origin mismatch - ",raster::origin(tmp)," ",firstOrigin)
          }
          if (!rgeos::gContainsProperly(mapshape, ei)) {
            temp <- system.time(
              tmp <- raster::mask(raster::crop(tmp, extent(mapshape),snap="near"),
                                  mapshape)
            )[[3]]          
            if (!silent) print(paste0("  boundary - masking time = ", temp))
          } else {
            if (!silent) print("  interior - not masked")
          }        
          r.list[[j]] <- tmp
          j <- j + 1
        } else {
          # no tile where expected
        }
      } else {
        print(paste0(  "lat=",lat," lon=",lon,"  exterior"))
      }
    }  
  }     
  if (!silent) print(warnings())
  if (j > 2) {
    if (!silent) print("calling merge")
    #m.sub <- do.call(merge, r.list))
    temp <- system.time(
      m.sub <- do.call(merge, c(r.list,list(tolerance=0.1)))
    )[3]
    if (!silent) print(temp)
  } else {
    m.sub <- r.list[[1]]
  }
  return(m.sub)
}
mapfilename <- function(lat,lon,resstr) {
  latChar <- ifelse(lat>=0,"n","s")
  lonChar <- ifelse(lon>=0,"e","w")
  return( paste0(latChar,sprintf("%02d",abs(lat)),"_",
                 lonChar,sprintf("%03d",abs(lon)),
                 resstr,".zip") )
}

