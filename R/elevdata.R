USStatevec <- NULL #c("PacificNorthWest","AK")
USParkvec <- NULL
CAProvincevec <- c("BC")
mapbuffer <- 2000     # meters
cropbox <- extent(-180, 0, -55, 60)

res3dplot <- 2500
loadStateElevs <-  FALSE
writeElevFile <- TRUE
drawRGL <- TRUE
saveRGL <- TRUE
drawPlotly <- TRUE
savePlotly <- TRUE
lon0to360=FALSE

datadir <- "c:/bda"                    #  trimmed elevation raster output data location
NAmericaDataDir <- "c:/bda/NorthAmerica"        #  zip file input subdirectories location
mapoutputdir <- "c:/bda/maps3d"        #  map file output location

#################################################################################
regionlist <- c("PacificNorthWest",
                "MountainWest",
                "DelMarVa",
                "NewEngland")
regionstates <- list(c("WA","OR","ID","MT"),
                     c("WA","OR","ID","MT","WY","CO","UT","NV"),
                     c("DE","MD","DC","VA"),
                     c("ME","NH","VT","MA","RI","CT")
                    )
#################################################################################
assign("last.warning", NULL, envir = baseenv())

library(tidyverse)
library(raster)
library(plotly)
library(htmlwidgets)
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(tigris)
library(sf)

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

addmapfiles <- function(filenames,lonChar,minLon,maxLon,latChar,minLat,maxLat) {
  #  first pass assume N and W quartersphere
  if (minLat < maxLat) {
    for (lat in seq(minLat,maxLat)) {
      if (minLon < maxLon) {
        for (lon in seq(minLon,maxLon))
          filenames <- c(filenames,
                         paste0(latChar,sprintf("%02d",lat),"_",lonChar,sprintf("%03d",lon),
                                "_1arc_v3_bil.zip"))
      }
    }
  }
  return(filenames)  
}


outputName <- paste0(c(USStatevec,CAProvincevec,USParkvec),collapse="-")
highelevation <- 3000

####    set up crop shape file
mapcrop <- NULL
statesInMap <- NULL
if (!is.null(USStatevec)) {
  USStatevec <- toupper(USStatevec)  # US State abbrev all upper case
  for (i in 1:length(regionlist)) {
    if (toupper(regionlist[[i]]) %in% USStatevec) {
      USStatevec <- setdiff(USStatevec,toupper(regionlist[[i]]))
      USStatevec <- union(USStatevec,toupper(regionstates[[i]]))
    }
  }
  statesInMap <- c(statesInMap,USStatevec)
  mcrop <- rgeos::gUnaryUnion(tigris::counties(statesInMap))
  if (mapbuffer > 0) 
    mcrop <- rgeos::gBuffer(sp::spTransform( mcrop, CRS( "+init=epsg:3857" ) ),
                            width=mapbuffer)
  mcrop <- sp::spTransform( mcrop, CRS( "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs" ) ) 
  if (is.null(mapcrop)) {
    mapcrop <- mcrop
  } else {
    mapcrop <- rgeos::gUnaryUnion(raster::union(mapcrop,mcrop))
  }
} 
if (!is.null(USParkvec)) {
  parkareas <- sf::st_read(paste0(datadir,
                                  "/ne_10m_parks_and_protected_lands",
                                  "/ne_10m_parks_and_protected_lands_area.shp"))
  
  if (is.null(mapcrop)) {
    mapcrop <- mcrop
  } else {
    mapcrop <- rgeos::gUnaryUnion(raster::union(mapcrop,mcrop))
  } 
}
if (!is.null(CAProvincevec)) {
  CAProvincevec <- unique(CAProvincevec)
  statesInMap <- c(statesInMap,CAProvincevec)
  canada <- raster::getData("GADM",country="CAN",level=1)
  mcrop <- rgeos::gSimplify(canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),],
                            tol=0.1)
  mcrop <- rgeos::gBuffer(sp::spTransform( mcrop, CRS( "+init=epsg:3857" ) ),
                          width=max(mapbuffer,16000))
  mcrop <- sp::spTransform( mcrop, CRS( "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs" ) ) 
  if (is.null(mapcrop)) {
    mapcrop <- mcrop
  } else {
    mapcrop <- rgeos::gUnaryUnion(raster::union(mapcrop,mcrop))
  } 
}
CP <- as(cropbox, "SpatialPolygons")
proj4string(CP) <- CRS(proj4string(mapcrop))
mapcrop <- gIntersection(mapcrop, CP, byid=TRUE)

plot(mapcrop)

j <- 1
r.list <- list()
if (loadStateElevs) {
  for (st in statesInMap) {
    load(paste0(datadir,"/",st,"elevs.rda"))
    r.list[[j]] <- elevations
    j <- j + 1
  }
  if (j > 2) {
    m.sub <- do.call(merge, r.list)
  } else {
    m.sub <- elevations
  }
} else {
  mapextent <- raster::extent(mapcrop)
  ELonMin <- floor(max(mapextent@xmin,0))
  ELonMax <- floor(max(mapextent@xmax,0))
  WLonMin <- ceiling(max(-mapextent@xmax,0))
  WLonMax <- ceiling(max(-mapextent@xmin,0))
  NLatMin <- floor(max(mapextent@ymin,0))
  NLatMax <- floor(max(mapextent@ymax,0))
  SLatMin <- ceiling(max(-mapextent@ymax,0))
  SLatMax <- ceiling(max(-mapextent@ymin,0))

  fn <- NULL
  if (NLatMin < NLatMax) {
    if (WLonMin < WLonMax) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"n",NLatMin,NLatMax)
    }
    if (ELonMin < ELonMax) {
      fn <- addmapfiles(fn,"e",WLonMin,WLonMax,"n",NLatMin,NLatMax)
    }
  }
  if (SLatMin < SLatMax) {
    if (WLonMin < WLonMax) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"s",SLatMin,SLatMax)
    }
    if (ELonMin < ELonMax) {
      fn <- addmapfiles(fn,"e",WLonMin,WLonMax,"s",SLatMin,SLatMax)
    }
  }
  tempd <- tempdir()
  rn <- gsub("_bil","",tools::file_path_sans_ext(fn))
  unfoundfn <- NULL
  for (i in 1:length(fn)) {
    if (file.exists(paste0(NAmericaDataDir,"/",fn[[i]]))) {
      print(paste0(NAmericaDataDir,"/",fn[[i]]))
      unzip(paste0(NAmericaDataDir,"/",fn[[i]]),exdir=tempd)
      tmp <- raster(paste0(tempd,"/",rn[[i]],".bil"))
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
      } else if (gIntersects(mapcrop, ei)) {
        print(paste0("boundary - masking time = ",system.time(
                 tmp <- raster::mask(raster::crop(tmp, extent(mapcrop)),
                                     mapcrop)
                          ))[[3]])
      } else {
        print("exterior - not used")
        tmp <- NULL
      }
      if (!is.null(tmp)) {
        r.list[[j]] <- tmp
        j <- j + 1
      }
    } else {
      print(paste0(NAmericaDataDir,"/",fn[[i]]," does not exist, ignored"))
      unfoundfn <- c(unfoundfn,fn[[i]])
    }
  }
  print(warnings())
  if (j > 2) {
    m.sub <- do.call(merge, r.list)
  } else {
    m.sub <- r.list[[1]]
  }
}
if (lon0to360) m.sub <- raster::rotate(m.sub)

if (writeElevFile & (length(statesInMap)==1)) {
  #  for now only write when a single state/province/whatever specified
  elevations <- raster::readAll(m.sub) 
  save(elevations,file=paste0(datadir,"/",statesInMap[[1]],"elevs.rda"))
} else {
  elevations <- m.sub
}

print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- max(1,floor(max(elevations@ncols,elevations@nrows)/res3dplot))
print(paste0("scaling raster down by a factor of ",sfact))
if (sfact > 1)
  elevations <- raster::aggregate(elevations,fact=sfact,fun=mean,
                           expand=TRUE,na.rm=TRUE)
print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))

if (drawRGL | drawPlotly) {

  yscale <- yRatio(elevations)
  elevations[is.na(elevations)] <- 0
  mmmelev <- raster::as.matrix(elevations)
  if (drawPlotly) {
    mmm <- mmmelev[,ncol(mmmelev):1]  #  flip east/west since row 1 is top
    ax <- list(title="longitude",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    ay <- list(title="latitude",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    az <- list(title="elevation",zeroline=FALSE,
               showline=FALSE,showticklabels=FALSE,showgrid=FALSE)
    p <- plotly::plot_ly(z = ~mmm,
                         colors = c("blue","yellow")) %>%
      plotly::add_surface(opacity=1.0) %>%
      plotly::layout(scene=list(xaxis=ax,yaxis=ay,zaxis=az,
                                aspectmode = "manual",
                                aspectratio = list(x=1,y=yscale,z=0.05),
                                camera=list(up=c(0,1,0),
                                            eye=c(0,1.25,0)) ) )
    p
    if (savePlotly) htmlwidgets::saveWidget(p,paste0(mapoutputdir,"/",outputName," map.html"))
  }
  if (drawRGL) {
    mmmelev <- raster::as.matrix(elevations)
    x <- seq(1,length.out=nrow(mmmelev))
    y <- seq(1,length.out=ncol(mmmelev))
    terrcolors <- colorRampPalette(c("blue","turquoise","aquamarine",
                                     "palegreen","yellowgreen",
                                     "chartreuse","greenyellow","green",
                                     "limegreen","forestgreen","darkgreen",
                                     "yellow","gold","goldenrod",
                                     "sienna","brown","gray75",
                                     "gray85","gray95","gray97","white"))(201)
    colidx <- floor(200*(mmmelev/highelevation)) + 1
    colidx[colidx>201] <- 201
    colidx[colidx<1] <- 1
    col <- terrcolors[colidx]

    rgl::par3d("windowRect"= c(100,100,1200,1000))
    userMatrix <- matrix(c(-0.02,-0.80,0.632,0,1,0,0.04,0,
                           -0.03,0.60,0.80,0,0,0,0,1),ncol=4,nrow=4)
    rgl::rgl.clear()
    rgl::surface3d(x,y,mmmelev,color=col)
    rgl::material3d(alpha=1.0,point_antialias=TRUE,smooth=TRUE,shininess=0)
    rgl::aspect3d(x=1,y=1/yscale,z=0.04)
    rgl::rgl.clear("lights")
    rgl::rgl.light(theta = 0, phi = 25,
                   viewpoint.rel=TRUE, specular="black")
    rgl::rgl.viewpoint(userMatrix=userMatrix,type="modelviewpoint")
    pan3d(2)  # right button for panning, doesn't play well with zoom)
    if (saveRGL) rgl::writeWebGL(dir=paste0(mapoutputdir), filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))
  }
}
