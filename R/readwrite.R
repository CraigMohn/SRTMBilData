loadSavedElevData <- function(savedNameVec,rasterDir,
                              noisy=FALSE,silent=FALSE) {
  j <- 1
  r.list <- list()
  for (st in savedNameVec) {
    fvec <- list.files(path=paste0(rasterDir,"/",st),
                       pattern=paste0(st,"elevs[0-9]{,2}.grd"))
    for (fn in fvec) {
      if (!silent) print(paste0("loading ",fn))
      elevations <- raster::raster(paste0(rasterDir,"/",st,"/",fn))
      if (noisy) print(elevations)
      r.list[[j]] <- elevations
      j <- j + 1
    }
  }
  if (!silent) print("calling merge")
  if (j > 2) {
    temp <- system.time(    
      m.sub <- do.call(merge, r.list)
    )[3]
    if (!silent) print(temp)
  } else {
    m.sub <- elevations
  }
  return(m.sub)
}
loadSavedFeatureData <- function(savedNameVec,rasterDir,
                                 noisy=FALSE,silent=FALSE) {
  lnames <- c("town","waterA","waterL","roads")
  rStack <- NULL
  for (layername in lnames) {
    j <- 1
    r.list <- list()
    for (st in savedNameVec) {
      fvec <- list.files(path=paste0(rasterDir,"/",st),
                         pattern=paste0(st,"features[0-9]{,2}_",layername,".grd"))
      for (fn in fvec) {
        if (!silent) print(paste0("loading ",fn))
        featureLayer <- raster::raster(paste0(rasterDir,"/",st,"/",fn))
        if (noisy) print(featureLayer)
        r.list[[j]] <- featureLayer
        j <- j + 1
      }
    }
    if (!silent) print("calling merge")
    if (j > 2) {
      temp <- system.time(    
        m.sub <- do.call(merge, r.list)
      )[3]
      if (!silent) print(temp)
    } else {
      m.sub <- featureLayer
    }
    if (is.null(rStack)) {
      rStack <- m.sub
    } else {
      rStack <- raster::addLayer(rStack,m.sub)      
    }
  }
  names(rStack) <- lnames
  return(rStack)
}
writeElevRaster <- function(elevations,maxrastercells,rasterDir,fname,
                            noisy=FALSE,silent=FALSE) {
  dir.create(file.path(rasterDir, fname))
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  if (!silent) print(paste0("saving raster data in ",nchunks," slices"))
  if (nchunks == 1) {
    writeRaster(elevations,file=paste0(rasterDir,"/",fname,"/",
                                       fname,"elevs.grd"),
                overwrite=TRUE)   
  } else {
    nrowchunk <- ceiling(raster::nrow(elevations)/nchunks)
    for (chunk in 1:nchunks) {
      if (noisy) print(paste0("cropping ",chunk))
      # overlap chunks by one row
      ylow <- raster::ymin(elevations) +  
              raster::yres(elevations)*((chunk-1)*nrowchunk)
      ylow[ylow < raster::ymin(elevations)] <- raster::ymin(elevations)
      yhi <-  raster::ymin(elevations) + 
              raster::yres(elevations)*(chunk*nrowchunk + 1)
      yhi[yhi > raster::ymax(elevations)] <- raster::ymax(elevations)
      chunkcrop <- raster::extent(raster::xmin(elevations),
                                  raster::xmax(elevations),
                                  ylow,yhi )
      chunkcrop <- as(chunkcrop, "SpatialPolygons")
      sp::proj4string(chunkcrop) <- sp::proj4string(elevations)
      chunkraster <- raster::trim(raster::crop(elevations,chunkcrop))
      if (noisy) print(paste0("writing ",chunk))
      writeRaster(chunkraster,
                  file=paste0(rasterDir,"/",fname,"/",
                              fname,"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  overwrite=TRUE)    
    }
  }
  return(NULL)
}
writeShapeFiles <- function(stname,shapefileDir,spTown,spRoads,spWaterA,spWaterL,
                            noisy=FALSE,silent=FALSE) {
  if (!silent) print(paste0("writing feature shapefiles for ",stname))
  raster::shapefile(spRoads,filename=paste0(shapefileDir,"/",stname,"Roads.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterA,filename=paste0(shapefileDir,"/",stname,"WaterA.shp"),
                    overwrite=TRUE)
  raster::shapefile(spWaterL,filename=paste0(shapefileDir,"/",stname,"WaterL.shp"),
                    overwrite=TRUE)
  raster::shapefile(spTown,filename=paste0(shapefileDir,"/",stname,"Town.shp"),
                    overwrite=TRUE)
  return(NULL)
}
writeFeatureRaster <- function(featureStack,maxrastercells,rasterDir,fname,
                               noisy=FALSE,silent=FALSE) {

  nchunks <- ceiling(raster::ncell(featureStack)/maxrastercells)
  rnrow <- raster::nrow(featureStack)
  nrowchunk <- ceiling(raster::nrow(featureStack)/nchunks)
  rxmin <- raster::xmin(featureStack)
  rxmax <- raster::xmax(featureStack)
  rymin <- raster::ymin(featureStack) 
  rymax <- raster::ymax(featureStack) 
  ryres <- raster::yres(featureStack)
  tempProj4 <- sp::proj4string(featureStack)
  if (!silent) print(paste0("saving featureraster data in ",nchunks," slices"))
  if (nchunks == 1) {
    if (noisy) print(paste0("writing ",fname,".grd"))
     writeRaster(featureStack,
                 file=paste0(rasterDir,"/",fname,"/",
                             fname,"features",
                             ".grd"),
                 bylayer=TRUE,suffix="names",   
                 datatype="INT1S",overwrite=TRUE)   
  } else {
    for (chunk in 1:nchunks) {
      ylow <- rymin + ryres*((chunk-1)*nrowchunk)
      ylow[ylow < rymin] <- rymin
      yhi <-  rymin + ryres*(chunk*nrowchunk + 1)
      yhi[yhi > rymax] <- rymax
      chunkcrop <- raster::extent(rxmin,rxmax,
                                  ylow,yhi)
      chunkcrop <- as(chunkcrop, "SpatialPolygons")
      sp::proj4string(chunkcrop) <- tempProj4
      if (noisy) print(paste0("cropping ",chunk))
      chunkraster <- raster::crop(featureStack,chunkcrop) 
      if (noisy) print(paste0("writing ",chunk))
      writeRaster(chunkraster,
                  file=paste0(rasterDir,"/",
                              fname,"features",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  bylayer=TRUE,suffix="names",   
                  datatype="INT1S",overwrite=TRUE)   
    }
  }
  return(NULL)
}
