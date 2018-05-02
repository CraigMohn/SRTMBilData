loadSavedElevData <- function(savedNameVec,rasterDir) {
  j <- 1
  r.list <- list()
  for (st in savedNameVec) {
    fvec <- list.files(path=rasterDir,
                       pattern=paste0(st,"elevs[0-9]{,2}.grd"))
    for (fn in fvec) {
      print(paste0("loading ",fn))
      elevations <- raster(paste0(datadir,"/rasterfiles/",fn))
      print(elevations)
      r.list[[j]] <- elevations
      j <- j + 1
    }
  }
  print("calling merge")
  if (j > 2) {
    m.sub <- do.call(merge, r.list)
  } else {
    m.sub <- elevations
  }
  return(m.sub)
}
loadSavedFeatureData <- function(savedNameVec,rasterDir) {
  lnames <- c("town","waterA","waterL","road")
  rStack <- NULL
  for (layername in lnames) {
    j <- 1
    r.list <- list()
    for (st in savedNameVec) {
      fvec <- list.files(path=rasterDir,
                         pattern=paste0(st,"features[0-9]{,2}_",layername,".grd"))
      for (fn in fvec) {
        print(paste0("loading ",fn))
        featureStack <- raster::stack(paste0(datadir,"/rasterfiles/",fn))
        print(featureStack)
        r.list[[j]] <- featureStack
        j <- j + 1
      }
    }
    print("calling merge")
    if (j > 2) {
      m.sub <- do.call(merge, r.list)
    } else {
      m.sub <- featureStack
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
writeElevRaster <- function(elevations,maxrastercells,rasterDir,fname) {
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  print(paste0("saving raster data in ",nchunks," slices"))
  if (nchunks == 1) {
    writeRaster(elevations,file=paste0(rasterDir,"/",
                                       fname,"elevs.grd"),
                overwrite=TRUE)   
  } else {
    nrowchunk <- ceiling(raster::nrow(elevations)/nchunks)
    for (chunk in 1:nchunks) {
      chunkcrop <- raster::extent(raster::xmin(elevations),
                                  raster::xmax(elevations),
                                  raster::ymin(elevations) + 
                                    raster::yres(elevations)*
                                    ((chunk-1)*nrowchunk - 10),
                                  raster::ymin(elevations) + 
                                    raster::yres(elevations)*
                                    (chunk*nrowchunk - 10)  )
      chunkcrop <- as(chunkcrop, "SpatialPolygons")
      sp::proj4string(chunkcrop) <- sp::proj4string(elevations)
      chunkraster <- raster::trim(raster::crop(elevations,chunkcrop))
      writeRaster(chunkraster,
                  file=paste0(rasterDir,"/",
                              fname,"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  overwrite=TRUE)    
    }
  }
  return(NULL)
}
writeFeatureRaster <- function(featureStack,maxrastercells,rasterDir,fname) {
  gc()  #  cleanup, this takes a lot of memory
  nchunks <- ceiling(raster::ncell(featureStack)/maxrastercells)
  print(paste0("saving raster data in ",nchunks," slices"))
  if (nchunks == 1) {
    writeRaster(featureStack,file=paste0(rasterDir,"/",
                                       fname,"features.grd"),
                bylayer=TRUE,suffix="names",   
                datatype="INT1S",overwrite=TRUE)   
  } else {
    nrowchunk <- ceiling(raster::nrow(featureStack)/nchunks)
    for (chunk in 1:nchunks) {
      chunkcrop <- raster::extent(raster::xmin(featureStack),
                                  raster::xmax(featureStack),
                                  raster::ymin(featureStack) + 
                                    raster::yres(featureStack)*
                                    ((chunk-1)*nrowchunk - 10),
                                  raster::ymin(featureStack) + 
                                    raster::yres(featureStack)*
                                    (chunk*nrowchunk - 10)  )
      chunkcrop <- as(chunkcrop, "SpatialPolygons")
      sp::proj4string(chunkcrop) <- sp::proj4string(featureStack)
      chunkraster <- raster::trim(raster::crop(featureStack,chunkcrop))
      writeRaster(chunkraster,
                  file=paste0(rasterDir,"/",
                              fname,"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  bylayer=TRUE,suffix="names",   
                  datatype="INT1S",overwrite=TRUE)   
    }
  }
  return(NULL)
}
writeShapeFiles <- function(stname,shapefileDir,spTown,spRoads,spWaterA,spWaterL) {
  print(paste0("writing feature shapefiles for ",stname))
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
readShapeFiles <- function(stname,shapefileDir) {
  print(paste0("loading features for ",stname))
  tmp <- raster::shapefile(paste0(shapefileDir,"/",stname,
                                  "Town.shp")) # SpatialPolygonsDF
  tmp <- sp::spTransform(tmp,workProj4)
  spTown <- tmp
  tmp <- raster::shapefile(paste0(shapefileDir,"/",stname,
                                  "Roads.shp")) # SpatialLinesDF
  tmp <- sp::spTransform(tmp,workProj4)
  spRoads <- tmp
  tmp <- raster::shapefile(paste0(shapefileDir,"/",stname,
                                  "WaterA")) # SpatialPolygonsDF
  tmp <- sp::spTransform(tmp,workProj4)
  spWaterA <- tmp
  tmp <- raster::shapefile(paste0(shapefileDir,"/",stname,
                                  "WaterL")) # SpatialLinesDF
  tmp <- sp::spTransform(tmp,workProj4)
  spWaterL <- tmp
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}

