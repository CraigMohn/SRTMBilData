#  rank change log
#  20180524 - add TIGER "S1500" (vehicle trails) to Roads class 2 
#  20180524 - add argument for TIGER data year, even state lines move 
#       so e sure sure years match when merging/stacking rasters

buildFeatureStack <- function(baseLayer,mapshape,
                              spList,
                              filterList=NULL,
                              maxRasterize=10000,
                              sliceFeatureBuffer=0,
                              polySimplify=0,polyMethod="vis",
                              polyWeighting=0.85,polySnapInt=0.0001  ) {
  #  baseLayer a rasterLayer
  #  returns a rasterStack  which has 4 layers which can be stored as ints
  
  #  collapse this down to spList and filterList
  bLrect <- as(raster::extent(baseLayer), "SpatialPolygons")
  sp::proj4string(bLrect) <- raster::crs(baseLayer)

  spList <- filterFeatures(spList,filterList)
  featureNames <- c("spTown","spRoads","spWaterA","spWaterL")
  for (x in featureNames) {
    if( x %in% names(spList) && 
        length(spList[[x]])>0 ) {
      assign(x,spList[[x]])
    } else {
      assign(x,NULL)
    }
  }
  if (length(spRoads)>0) {
    print("roads")
    spRoads <- sxdfMask(spRoads,bLrect)
  }
  if (length(spRoads)>0) {
    print(paste0(nrow(spRoads)," roads to process"))
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
  if (length(spWaterL)>0) {
    print("water lines")
    spWaterL <- sxdfMask(spWaterL,bLrect)
  }
  if (length(spWaterL)>0) {
    print(paste0(nrow(spWaterL)," water lines to process"))
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
  if (length(spWaterA)>0) {
    print("water polygons")
    if (sliceFeatureBuffer>0) {
      bLrectA <- bufferUnion(bLrect,
                              mapbuffer=sliceFeatureBuffer,
                              mapunion=NULL)
    } else {
      bLrectA <- bLrect
    }
    spWaterA <- sxdfMask(spWaterA,bLrectA)
  }
  if (length(spWaterA)>0) {
    print(paste0(nrow(spWaterA)," water areas to process"))
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
  if (length(spTown)>0) {
    print("towns")
    if (sliceFeatureBuffer>0) {
      bLrectA <- bufferUnion(bLrect,
                             mapbuffer=sliceFeatureBuffer,
                             mapunion=NULL)
    } else {
      bLrectA <- bLrect
    }
    spTown <- sxdfMask(spTown, bLrectA)
  }
  if (length(spTown)>0) {
    print(paste0(nrow(spTown)," towns to process"))
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
  s <- raster::stack(tlayer,wAlayer,wLlayer,rlayer)
  names(s) <- c("town","waterA","waterL","roads")  
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
  sp::proj4string(CP) <- raster::crs(templateRaster)
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
      print(paste0("  ",round(system.time(
        rlayer$rasterize(spdf=sxdf[first:last,],field="value",background=0)
      )[[3]],digits=2)," rasterizing"))
      if (nloops>1) {
        print(paste0("  ",round(system.time(
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

