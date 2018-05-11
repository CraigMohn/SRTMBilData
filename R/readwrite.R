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
    print(system.time(    
      m.sub <- do.call(merge, r.list)
    )[3])
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
        featureLayer <- raster::raster(paste0(datadir,"/rasterfiles/",fn))
        print(featureLayer)
        r.list[[j]] <- featureLayer
        j <- j + 1
      }
    }
    print("calling merge")
    if (j > 2) {
      m.sub <- do.call(merge, r.list)
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
      print(paste0("cropping ",chunk))
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
      print(paste0("writing ",chunk))
      writeRaster(chunkraster,
                  file=paste0(rasterDir,"/",
                              fname,"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  overwrite=TRUE)    
    }
  }
  return(NULL)
}
writeFeatureRaster <- function(featureStack=NULL,elevations=NULL,mapshape=NULL,
                               maxrastercells,rasterDir,fname,
                               spTown=NULL,spRoads=NULL,
                               spWaterA=NULL,spWaterL=NULL,
                               polySimplify=0.0,polyMethod="vis", 
                               polyWeighting=0.85,polySnapInt=0.0001) {
  gc()        #  cleanup, this takes a lot of memory
  
  if (is.null(featureStack) & (is.null(elevations)|is.null(mapshape)))
    stop("no feature rasterStack and no elevation/mapshape pair")
  if (!is.null(featureStack)) {
    nchunks <- ceiling(raster::ncell(featureStack)/maxrastercells)
    rnrow <- raster::nrow(featureStack)
    nrowchunk <- ceiling(raster::nrow(featureStack)/nchunks)
    rxmin <- raster::xmin(featureStack)
    rxmax <- raster::xmax(featureStack)
    rymin <- raster::ymin(featureStack) 
    rymax <- raster::ymax(featureStack) 
    ryres <- raster::yres(featureStack)
    tempProj4 <- sp::proj4string(featureStack)
  } else {
    nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
    rnrow <- raster::nrow(elevations)
    nrowchunk <- ceiling(raster::nrow(elevations)/nchunks)
    rxmin <- raster::xmin(elevations)
    rxmax <- raster::xmax(elevations)
    rymin <- raster::ymin(elevations) 
    rymax <- raster::ymax(elevations) 
    ryres <- raster::yres(elevations)
    tempProj4 <- sp::proj4string(elevations)
  }
  print(paste0("saving featureraster data in ",nchunks," slices"))
  if (nchunks == 1) {
    if (is.null(featureStack))
      featureStack <- buildFeatureStack(elevations,mapshape=mapshape,
                                        spTown=spTown,spWaterA=spWaterA,
                                        spWaterL=spWaterL,spRoads=spRoads,
                                        polySimplify=polySimplify,
                                        polyMethod=polyMethod, 
                                        polyWeighting=polyWeighting,
                                        polySnapInt=polySnapInt) %>%
            raster::trim(.)
        writeRaster(featureStack,file=paste0(rasterDir,"/",
                                       fname,"features.grd"),
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
      if (is.null(featureStack)) {
        print(paste0("cropping, building and masking ",chunk))
        chunkelev <- raster::crop(elevations,chunkcrop)
        chunkraster <- buildFeatureStack(chunkelev,mapshape=mapshape,
                                         spTown=spTown,spWaterA=spWaterA,
                                         spWaterL=spWaterL,spRoads=spRoads,
                                         polySimplify=polySimplify,
                                         polyMethod=polyMethod, 
                                         polyWeighting=polyWeighting,
                                         polySnapInt=polySnapInt) %>%
                  raster::trim(.)     
      } else {
        print(paste0("cropping ",chunk))
        chunkraster <- raster::crop(featureStack,chunkcrop) %>%
                  raster::trim(.)              
      }
      print(paste0("writing ",chunk))
      print(chunkraster)
      plot(chunkraster[["road"]])
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
readShapeFiles <- function(stname,shapefileDir,workProj4) {
  print(paste0("loading shapefiles for ",stname))
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
  print("done loading shapefiles")
  return(list(spTown=spTown,spRoads=spRoads,spWaterA=spWaterA,spWaterL=spWaterL))
}

