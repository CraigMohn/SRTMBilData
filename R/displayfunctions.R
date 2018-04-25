draw3DMapRgl <- function(elevations,mapPolygon,
                         spTown,spRoads,spWaterA,spWaterL,
                         waterLevel,colors="default",
                         maxElev=3000,vScale=1.5,   
                         rglNAcolor="Blue",citycolor="White",rglNegcolor=NA,
                         saveRGL=FALSE,mapoutputdir=NA,outputName=NA,
                         fastAreas=TRUE,polySimplify=0,
                         polyMethod=NULL,polyWeighting=0.85,
                         polySnapInt=0.0001) {
  
  yscale <- yRatio(elevations)
  mmmelev <- raster::as.matrix(elevations)
  x <- seq(1,length.out=nrow(mmmelev))
  y <- seq(1,length.out=ncol(mmmelev))

  CP <- as(doubleExtent(elevations), "SpatialPolygons")
  sp::proj4string(CP) <- CRS(sp::proj4string(mapPolygon))
  
  
  if (!(is.null(spTown) &
        is.null(spRoads) &
        is.null(spWaterA) &
        is.null(spWaterL))) {    
    
    print("building features")
    featurelayer <- raster::raster(nrow=nrow(elevations),ncol=ncol(elevations), 
                                   ext=extent(elevations),crs=crs(elevations))
    featurelayer[] <- 0
    if (!is.null(spTown)) {
      print("towns")
      tspTown <- raster::intersect(spTown, mapPolygon)
      plot(tspTown,col=citycolor)
      if (!is.null(tspTown)) {
        if (fastAreas) {
          featurelayer <- raster::rasterize(tspTown,featurelayer,
                                            field=1,update=TRUE)
          
        } else {
          for (i in 1:nrow(tspTown)) {
            cat(paste0("\r",i,"  ",tspTown@data[i,"NAME"],
                       "                                    "))
            #plot(tspTown[i,],col=citycolor)
            featurelayer <- raster::rasterize(tspTown[i,],featurelayer,
                                              field=1,update=TRUE)
          }
        }
        if (!is.finite(featurelayer@data@min)) {
          warning("featurelayer mess-up")
          print(featurelayer@data@min)
        }
        print(paste0(nrow(tspTown)," towns drawn"))
      } else {
        print("no towns to add")
      }
    }
    if (!is.null(spWaterA)) {
      print("water polygons")
      tspWaterA <- filterWaterA(spWaterA,level=waterLevel) 
      if (!is.null(tspWaterA)) 
        tspWaterA <- raster::intersect(tspWaterA, mapPolygon)
      if (!is.null(tspWaterA)) {
        if (fastAreas) {
          nWaterA <- nrow(tspWaterA)
          featurelayer <- raster::rasterize(tspWaterA,featurelayer,
                                            field=2,update=TRUE)
          
        } else {
          nWaterA <- 0
          print(paste0(nrow(tspWaterA)," water areas to process"))
          for (i in 1:nrow(tspWaterA)) {
            if (polySimplify == 0) {
              temp <- rgeos::gSimplify(tspWaterA[i,],tol=0,topologyPreserve=TRUE)
            } else {
              temp <- rmapshaper::ms_simplify(tspWaterA[i,],keep=polySimplify,
                                              method=polyMethod,
                                              weighting=polyWeighting,
                                              snap=TRUE,snap_interval=polySnapInt)
              temp <- sp::spTransform(temp,sp::CRS(sp::proj4string(mapPolygon)))
            }
             if (rgeos::gIntersects(temp,mapPolygon)) {
              cat(paste0("\r",i,"  ",tspWaterA@data[i,"NAME"],
                         "                                 "))
              temp <- rgeos::gIntersection(temp,mapPolygon,drop_lower_td=TRUE)
              featurelayer <- raster::rasterize(temp,featurelayer,
                                                field=2,update=TRUE)
              nWaterA <- nWaterA + 1
            }
          }
        }
        plot(tspWaterA, col="blue")
        if (!is.finite(featurelayer@data@min)) {
          warning("featurelayer mess-up")
          print(featurelayer@data@min)
        }
        print(paste0(nWaterA," water polygons drawn"))
      } else {
        print("no water polygons to add")
      }
    }
    if (!is.null(spWaterL)) {
      print("water lines")
      tspWaterL <- filterWaterL(spWaterL,level=waterLevel)
      if (!is.null(tspWaterL))tspWaterL <- raster::intersect(tspWaterL, mapPolygon)
      if (!is.null(tspWaterL)) {
        nWaterL <- nrow(tspWaterL)
        tspWaterL <- rgeos::gLineMerge(tspWaterL)
        plot(tspWaterL, col="blue")
        featurelayer <- raster::rasterize(tspWaterL,featurelayer,
                                          field=2,update=TRUE)
        if (!is.finite(featurelayer@data@min)) {
          warning("featurelayer mess-up")
          print(featurelayer@data@min)
        }
        print(paste0(nWaterL," water lines drawn"))
      } else {
        print("no water lines to add")
      }
    }
    if (!is.null(spRoads)) {
      print("roads")
      tspRoads <- raster::intersect(spRoads, mapPolygon) 
      if (!is.null(tspRoads)) {
        nRoads <- nrow(tspRoads)
        tspRoads <- rgeos::gLineMerge(tspRoads)
        plot(tspRoads)
        featurelayer <- raster::rasterize(tspRoads,featurelayer,
                                          field=3,update=TRUE)
        if (!is.finite(featurelayer@data@min)) {
          warning("featurelayer mess-up")
          print(featurelayer@data@min)
        }
        print(paste0(nRoads," roads drawn"))
      } else {
        print("no roads to add")
      }
    }
    features <- raster::as.matrix(featurelayer)
    print("features done")
  }
  if (colors == "beach") {
    terrcolors <- colorRampPalette(c("blue","bisque1","bisque2","bisque3",
                                     "palegreen","greenyellow","lawngreen",
                                     "chartreuse","green","springgreen",
                                     "limegreen","forestgreen","darkgreen",
                                     "olivedrab","darkkhaki","darkgoldenrod",
                                     "sienna","brown","saddlebrown","rosybrown",
                                     "gray35","gray45","gray55",
                                     "gray65","gray70","gray75","gray85"))(206)
  } else {
    terrcolors <- colorRampPalette(c("blue","darkturquoise","turquoise","aquamarine",
                                     "palegreen","greenyellow","lawngreen",
                                     "chartreuse","green","springgreen",
                                     "limegreen","forestgreen","darkgreen",
                                     "olivedrab","darkkhaki","darkgoldenrod",
                                     "sienna","brown","saddlebrown","rosybrown",
                                     "gray35","gray45","gray55",
                                     "gray65","gray70","gray75","gray85"))(206)
  }
  plot(rep(1,206),col=terrcolors, pch=19,cex=2)
  
  tmpelev <- mmmelev/maxElev # don't worry about memory, not the constraint here
  tmpelev[is.na(tmpelev)] <- 0
  tmpelev <- sign(tmpelev)*sqrt(abs(tmpelev)) # f(0)=0, f(1)=1, f'(x>0) decreasing, reasonable for x<0
  tmpelev[mmmelev <= 0] <- 0  #  go off original 
  colidx <- floor(200*tmpelev) + 1
  colidx[colidx>201] <- 201
  colidx[colidx<1] <- 1
  colidx <- colidx + 5
  colidx[mmmelev == 0]  <- 1
  col <- terrcolors[colidx]
  if (!is.na(rglNegcolor)) col[mmmelev < -10] <- gplots::col2hex(rglNegcolor)
  col[is.na(mmmelev)] <- gplots::col2hex(rglNAcolor)
  mmmelev[is.na(mmmelev)] <- -50
  if (showCities | showWater | showRoads) {
    col[features==1] <- gplots::col2hex(citycolor)
    col[features==2] <- gplots::col2hex("Blue")
    col[features==3] <- gplots::col2hex("Black")
  }
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
  if (saveRGL) 
    rgl::writeWebGL(dir=paste0(mapoutputdir), 
                    filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))

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
