
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


sxdfMask <- function(sxdf,poly,keepTouch=FALSE) {
  
  #  return NULL if a) either is NULL or b) no overlap
  #  horrible multiple returns, but....
  if (is.null(sxdf) | is.null(poly)) return(NULL) 
  
  tmpdata <- sxdf@data
  tmp.1 <- rgeos::gIntersects(sxdf, poly, byid=TRUE)
  tmp.2 <- as.logical(apply(tmp.1, 2, function(x) {sum(x)} ))
  if (sum(tmp.2) == 0) return(NULL)
  
  #  keep only intersecting sp objects in dataframe
  tmpgeo <- sxdf[tmp.2,]
  tmpdata <- tmpdata[tmp.2,]
  row.names(tmpgeo) <- row.names(tmpdata)
  if (!keepTouch) {
    return(raster::intersect(tmpgeo,poly))
  } else if (class(sxdf)=="SpatialLinesDataFrame") {
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

