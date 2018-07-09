bufferUnion <- function(spObj,mapbuffer,mapunion,
                        outCRS="",
                        bufferCRS="+init=epsg:3857",
                        capStyle="FLAT",
                        joinStyle="BEVEL",
                        simplifytol=0) {
  if (missing(outCRS)){
    outCRS <- sp::proj4string(spObj)
  }
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
  } else if (is.null(b)) {
    return(a)
  } else {
    return(rbind(a,b))
  }
}
rbindList <- function(spListA,spListB) {
  out <- list()
  for (x in union(names(spListA),names(spListB))) {
    out[[x]] <- rbind_NULLok(spListA[[x]],spListB[[x]])
  }
  return(out)
}
doubleExtent <- function(objectWithExtent) {
  bb <- sp::bbox(objectWithExtent)
  return(extent(max(-180,(3*bb[1,1]-bb[1,2])/2),
                min(180,(-bb[1,1]+3*bb[1,2])/2),
                max(-90,(3*bb[2,1]-bb[2,2])/2),
                min(90,(-bb[2,1]+3*bb[2,2])/2)))
}
quickmask <- function(r,sp,rectangle=FALSE) {
  if (rectangle) {
    return(raster::crop(r,sp))
  } else {
    return(raster::mask(r,sp))
  }
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
argCaseFix <- function(argchar,validCase) {
  for (x in validCase) {
    if (toupper(argchar)==toupper(x)) return(x)
  }
}
UTMProj4 <- function(lon,lat) {
  hemisphere <- ifelse(lat > 0,"+north","+south")
  return(paste0("+proj=utm +zone=",UTMzone(lon)," ",hemisphere,
         " +ellps=WGS84 +datum=WGS84"))
}
UTMzone <- function(lon) {
  #  from stackoverflow answer by Josh O'Brien
  return((floor((lon + 180)/6) %% 60) + 1) 
}
CRS_LambertAzimuthalEqualArea <- function() {
  "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs "
}
CRS_AlbersEqualArea <- function() {
  "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
}

