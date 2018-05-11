draw3DMapTrack <- function(mapRaster,trackdf=NULL,
                           featureLevels=c(3,3,2,5), #towns,roads,waterA,waterL
                           maxElev=3000,vScale=1.5,
                           colors="default",
                           citycolor="White",roadcolor="Black",
                           watercolor="Blue",glaciercolor="White",
                           rglNAcolor="Blue",rglNegcolor=NA,
                           trackcolor="Magenta",trackCurve=FALSE,
                           trackCurveElevFromRaster=TRUE,trackCurveHeight=10,
                           saveRGL=FALSE,mapoutputdir=NA,
                           outputName="most recent") {

  if (class(mapRaster)=="RasterLayer") {
    elevations <- mapRaster
  } else {
    elevations <- mapRaster[["elevations"]]
  }  
  mmmelev <- raster::as.matrix(elevations)
  x <- seq(1,length.out=nrow(mmmelev))
  y <- seq(1,length.out=ncol(mmmelev))
  yscale <- yRatio(elevations)

  terrcolors <- terrainColors(colors,206)
                              
  #  assign elevation-based colors 
  tmpelev <- mmmelev/maxElev    #  rescale in terms of maximum
  tmpelev[is.na(tmpelev)] <- 0
  #tmpelev <- sign(tmpelev)*sqrt(abs(tmpelev)) # f(0)=0, f(1)=1, f'(x>0) decreasing, reasonable for x<0
  tmpelev[mmmelev <= 0] <- 0  #  go off original 
  colidx <- floor(200*tmpelev) + 1
  colidx[colidx>201] <- 201   #  cap at maxElev
  colidx[colidx<1] <- 1       #  and at 0
  colidx <- colidx + 5
  colidx[mmmelev == 0]  <- 1
  col <- terrcolors[colidx]
  if (!is.na(rglNegcolor)) col[mmmelev < -10] <- gplots::col2hex(rglNegcolor)
  col[is.na(mmmelev)] <- gplots::col2hex(rglNAcolor)
  mmmelev[is.na(mmmelev)] <- -20    #  have missing elevations slightly below zero
  
  #  draw cities, water and roads in that order
  if ("town" %in% names(mapRaster)) {
    town <- as.matrix(mapRaster[["town"]])
    col[ town >= featureLevels[1] ] <- gplots::col2hex(citycolor)
    town <- NULL
  }
  if ("waterA" %in% names(mapRaster)) {
    waterA <- as.matrix(mapRaster[["waterA"]])
    col[ waterA >= featureLevels[3] ] <- gplots::col2hex(watercolor)
  }
  if ("waterL" %in% names(mapRaster)) {
    waterL <- as.matrix(mapRaster[["waterL"]])
    col[ waterL >= featureLevels[4] ] <- gplots::col2hex(watercolor)
    waterL <- NULL    
  }
  if ("waterA" %in% names(mapRaster)) {
      col[ waterA == 8 ] <- gplots::col2hex(glaciercolor) # Glaciers overwrite other water
      waterA <- NULL
  }
  if ("road" %in% names(mapRaster)) {
    road <- as.matrix(mapRaster[["road"]])
    col[ road >= featureLevels[2] ] <- gplots::col2hex(roadcolor)
    road <- NULL
  }
  if (!is.null(trackdf)) {
    if (!trackCurve) {
      cellnum <- cellFromXY(mapRaster,cbind(trackdf$lon,trackdf$lat))
      col[cellnum] <- trackcolor
    } else {
      xmin <- raster::extent(mapRaster)[1]
      xmax <- raster::extent(mapRaster)[2]
      xmin <- raster::extent(mapRaster)[3]
      xmax <- raster::extent(mapRaster)[4]
      xlen <-
        (raster::pointDistance(cbind(xmin,ymin),cbind(xmax,ymin),lonlat=TRUE) +
           raster::pointDistance(cbind(xmin,ymax),cbind(xmax,ymax),lonlat=TRUE)) /
        2
      ylen <-
        (raster::pointDistance(cbind(xmin,ymin),cbind(xmin,ymax),lonlat=TRUE) +
           raster::pointDistance(cbind(xmax,ymin),cbind(xmax,ymax),lonlat=TRUE)) /
        2
      
      xpath <- xlen * (1 - (trackdf$lat-ymin)/(ymax-ymin))
      ypath <- ylen * (trackdf$lon-xmin)/(xmax-xmin)
      if (trackCurveElevFromRaster)  {
        zpath <- raster::extract(mapRaster,cbind(trackdf$lon,trackdf$lat),
                             method="bilinear") +  trackCurveHeight
      } else {
        zpath <- trackdf$altitude.m  +  trackCurveHeight
      }
    }
  }
  
  #  and output the graph using rgl
  rgl::par3d("windowRect"= c(100,100,1200,1000))
  userMatrix <- matrix(c(-0.02,-0.80,0.632,0,1,0,0.04,0,
                         -0.03,0.60,0.80,0,0,0,0,1),ncol=4,nrow=4)
  rgl::rgl.clear()
  rgl::surface3d(x,y,mmmelev,color=col)
  rgl::material3d(alpha=1.0,point_antialias=TRUE,smooth=TRUE,shininess=0)
  rgl::aspect3d(x=1,y=1/yscale,z=0.035*vScale)
  rgl::rgl.clear("lights")
  rgl::rgl.light(theta = 0, phi = 15,
                 viewpoint.rel=TRUE, specular="black")
  rgl::rgl.viewpoint(userMatrix=userMatrix,type="modelviewpoint")
  pan3d(2)  # right button for panning, doesn't play well with zoom)
  if (!is.null(trackdf) & trackCurve) {
    rgl::lines3d(xpath,ypath,zpath,
                 lwd=2,col = trackcolor, alpha=1.0)
    
  }
  if (saveRGL) 
    rgl::writeWebGL(dir=paste0(mapoutputdir), 
                    filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))
  return(NULL)
}
terrainColors <- function(palettename="default",numshades=206) {
  if (palettename == "default") {
    terrcolors <- 
      colorRampPalette(c("blue","darkturquoise","turquoise","aquamarine",
                         "palegreen","greenyellow","lawngreen",
                         "chartreuse","green","springgreen",
                         "limegreen","forestgreen","darkgreen",
                         "olivedrab","darkkhaki","darkgoldenrod",
                         "sienna","brown","saddlebrown","rosybrown",
                         "gray35","gray45","gray55",
                         "gray65","gray70","gray75","gray85"))(numshades)
    
  } else if (palettename == "beach") {
    terrcolors <- 
      colorRampPalette(c("blue","bisque1","bisque2","bisque3",
                         "palegreen","greenyellow","lawngreen",
                         "chartreuse","green","springgreen",
                         "limegreen","forestgreen","darkgreen",
                         "olivedrab","darkkhaki","darkgoldenrod",
                         "sienna","brown","saddlebrown","rosybrown",
                         "gray35","gray45","gray55",
                         "gray65","gray70","gray75","gray85"))(numshades)
  } else if (palettename == "viridis") {
    terrcolors <- 
      viridis::viridis_pal(begin=0.2,end=0.9,direction=1,option="D")(numshades)
  } else if (palettename == "plasma") {
    terrcolors <- 
      viridis::viridis_pal(begin=0.0,end=1.0,direction=-1,option="D")(numshades)
  }
  return(terrcolors)
}
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
