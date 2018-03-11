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

mapLibSelector <- 1
mapWindow <- NULL # <- c(-122.64,-121.67,45.95,46.70) # NULL
USStatevec <- NULL # <- c("MountainWest","CA","NM","AZ")
USParkvec <- c("CRLA") 
CAProvincevec <- NULL #<- c("BC","AB","SK")
worldCountryvec <- NULL # <- "AUS" #NULL # <- c("ESP","PRT","FRA") # http://kirste.userpage.fu-berlin.de/diverse/doc/ISO_3166.html
mapbuffer <- 3000     # meters, expands overall area
mapmergebuffer <- 20  # meters, expands categories before merging with others
#  southern slice of BC,AB,SK 
cropbox <- raster::extent(-180, 170, -50, 52.1)
res3dplot <- 5000
loadStateElevs <-  FALSE
writeElevFile <- FALSE
forceRes <- NULL
maxrastercells <- 500000000
highelevation <- 3000
vertscale <- 1.0
drawRGL <- TRUE
saveRGL <- TRUE
drawPlotly <- FALSE
savePlotly <- FALSE
lon0to360=FALSE

datadir <- "c:/bda"                    #  base dir - subs include shapefiles, rasterfile 
mapoutputdir <- "c:/bda/maps3d"        #  map file output location
NAmericaDataDir <- "c:/bda/NorthAmerica"   #  zip file input subdirectories location
EuropeDataDir <- "c:/bda/Europe 3s"        #  zip file input subdirectories location
AustraliaDataDir <- "c:/bda/Australia 3s"        #  zip file input subdirectories location
mapDataDir <- c(NAmericaDataDir,EuropeDataDir,AustraliaDataDir)[[mapLibSelector]]
resstr <- c("_1arc_v3_bil","_3arc_v2_bil","_3arc_v2_bil")[[mapLibSelector]]


#################################################################################
regionlist <- c("PacificNorthWest",
                "MountainWest",
                "DelMarVa",
                "NewEngland",
                "Maritimes")
regionstates <- list(c("WA","OR","ID","MT"),
                     c("WA","OR","ID","MT","WY","CO","UT","NV"),
                     c("DE","MD","DC","VA"),
                     c("ME","NH","VT","MA","RI","CT"),
                     c("NS","PE","NB","NL")
                    )
#################################################################################
assign("last.warning", NULL, envir = baseenv())


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

addmapfiles <- function(filenames,lonChar,minLon,maxLon,latChar,minLat,maxLat,
                        resstr="_1arc_v3_bil") {
  #  first pass assume N and W quartersphere
  if (minLat <= maxLat) {
    for (lat in seq(minLat,maxLat)) {
      if (minLon <= maxLon) {
        for (lon in seq(minLon,maxLon))
          filenames <- c(filenames,
                         paste0(latChar,sprintf("%02d",lat),"_",lonChar,sprintf("%03d",lon),
                                resstr,".zip"))
      }
    }
  }
  return(filenames)  
}
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


outputName <- paste0(c(USStatevec,CAProvincevec,USParkvec),collapse="-")

####    set up crop shape file
mapcrop <- NULL
statesInMap <- NULL
if (!is.null(mapWindow)) {
  mapcrop <- raster::extent(mapWindow)
  CP <- as(mapcrop, "SpatialPolygons")
  sp::proj4string(CP) <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"
  mapcrop <- rgeos::gUnaryUnion(CP)
}
if (!is.null(USStatevec)) {
  USStatevec <- toupper(USStatevec)  # US State abbrev all upper case
  for (i in 1:length(regionlist)) {
    if (toupper(regionlist[[i]]) %in% USStatevec) {
      USStatevec <- setdiff(USStatevec,toupper(regionlist[[i]]))
      USStatevec <- union(USStatevec,toupper(regionstates[[i]]))
    }
  }
  statesInMap <- union(statesInMap,USStatevec)
  mcrop <- rgeos::gUnaryUnion(tigris::counties(statesInMap)) #sf dataframe
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
}
if (!is.null(USParkvec)) {
  USParkvec <- unique(toupper(USParkvec))
  statesInMap <- union(statesInMap,USParkvec)
  pdir <- "nps_boundary"
  pfile <- paste0(pdir,".shp")
  parkareas <- sf::st_read(paste0(datadir,"/",pdir,"/",pfile))  # sf dataframe
  #parkareas <- rgeos::gUnaryUnion(parkareas)
  #  parknames <- parkareas[,c("UNIT_NAME","UNIT_CODE")]
  #  sf::st_geometry(parknames) <- NULL
  #  parknames <- parknames[order(parknames$UNIT_CODE),]
  #  write.csv(parknames,paste0(datadir,"/parknames.csv"))
  parkareas <- parkareas[parkareas$UNIT_CODE %in% USParkvec,"geometry"]
  mcrop <- as(parkareas, "Spatial")  
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
}

if (!is.null(CAProvincevec)) {
  CAProvincevec <- unique(toupper(CAProvincevec))
  for (i in 1:length(regionlist)) {
    if (toupper(regionlist[[i]]) %in% CAProvincevec) {
      CAProvincevec <- setdiff(CAProvincevec,toupper(regionlist[[i]]))
      CAProvincevec <- union(CAProvincevec,toupper(regionstates[[i]]))
    }
  }
  statesInMap <- union(statesInMap,CAProvincevec)
  canada <- raster::getData("GADM",country="CAN",level=1) # raster + spatial
  mcrop <- rgeos::gUnaryUnion(canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),])
  # simplify - BC coast is extremely complex
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1)
}
if (!is.null(worldCountryvec)) {
  worldCountryvec <- unique(toupper(worldCountryvec))
  statesInMap <- union(statesInMap,worldCountryvec)
  mcrop <- NULL
  for (c in worldCountryvec) {
    cmap <- raster::getData("GADM",country=worldCountryvec,level=0) # raster + spatial
    cmap <- rgeos::gUnaryUnion(cmap)
    if (is.null(mcrop)) {
      mcrop <- cmap
    } else {
      mcrop<- rgeos::gUnaryUnion(raster::union(mcrop,cmap))
    }
  }
  mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1)
}
mapcrop <- bufferUnion(mapcrop,mapbuffer=mapbuffer,NULL,simplifytol = 0)
CP <- as(cropbox, "SpatialPolygons")
sp::proj4string(CP) <- CRS(sp::proj4string(mapcrop))
mapcrop <- rgeos::gIntersection(mapcrop, CP, byid=TRUE)

plot(mapcrop)

j <- 1
r.list <- list()
if (loadStateElevs) {
  for (st in statesInMap) {
    fvec <- list.files(path=paste0(datadir,"/rasterfiles"),
                        pattern=paste0(st,"elevs[0-9]{,2}.grd"))
    for (fn in fvec) {
      print(paste0("loading ",fn))
      elevations <- raster(paste0(datadir,"/rasterfiles/",fn))
      r.list[[j]] <- elevations
      j <- j + 1
    }
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
  if (NLatMax > 0) {
    if (WLonMax > 0) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"n",NLatMin,NLatMax,resstr)
    }
    if (ELonMax > 0) {
      fn <- addmapfiles(fn,"e",ELonMin,ELonMax,"n",NLatMin,NLatMax,resstr)
    }
  }
  if (SLatMax > 0) {
    if (WLonMax > 0) {
      fn <- addmapfiles(fn,"w",WLonMin,WLonMax,"s",SLatMin,SLatMax,resstr)
    }
    if (ELonMax > 0) {
      fn <- addmapfiles(fn,"e",ELonMin,ELonMax,"s",SLatMin,SLatMax,resstr)
    }
  }
  tempd <- tempdir()
  rn <- gsub("_bil","",tools::file_path_sans_ext(fn))
  unfoundfn <- NULL
  firstRes <- NULL
  firstProj <- NULL
  firstOrigin <- NULL
  for (i in 1:length(fn)) {
    if (file.exists(paste0(mapDataDir,"/",fn[[i]]))) {
      cat("\n")
      print(paste0(mapDataDir,"/",fn[[i]]))
      unzip(paste0(mapDataDir,"/",fn[[i]]),exdir=tempd)
      tmp <- raster(paste0(tempd,"/",rn[[i]],".bil"))
      if (is.null(firstOrigin)) {
        print(tmp)
        firstOrigin <- raster::origin(tmp)
        print(firstOrigin)
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
      if (!identical(firstRes,raster::res(tmp))) {
        print("resolution differs from the first tile - resampling")
        print(raster::res(tmp))
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
        #raster::res(newraster) <- firstRes
        raster::origin(newraster) <- firstOrigin
        print(tmp)  
        tmp <- raster::resample(tmp, newraster)
        print(tmp)
      } else {
        if (max(abs(firstOrigin-raster::origin(tmp))) > 0.0000001)  
          warning("origin mismatch - ",raster::origin(tmp)," ",firstOrigin)
      }
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
                 tmp <- raster::mask(raster::crop(tmp, extent(mapcrop),snap="near"),
                                     mapcrop)
                          ))[[3]])
        print(origin(tmp))
      } else {
        print("exterior - not used")
        tmp <- NULL
      }
      if (!is.null(tmp)) {
        r.list[[j]] <- tmp
        j <- j + 1
      }
    } else {
      print(paste0(mapDataDir,"/",fn[[i]]," does not exist, ignored"))
      unfoundfn <- c(unfoundfn,fn[[i]])
    }
  }
  print(warnings())
  print("calling merge")
  if (j > 2) {
    #m.sub <- do.call(merge, r.list))
    m.sub <- do.call(merge, c(r.list,list(tolerance=0.1)))
  } else {
    m.sub <- r.list[[1]]
  }
}
if (lon0to360) m.sub <- raster::rotate(m.sub)
elevations <- m.sub

if (writeElevFile & (length(statesInMap)==1)) {
  #  for now only write when a single state/province/whatever specified
  #writeRaster(elevations,file=paste0(datadir,"/",statesInMap[[1]],"elevs.grd"),
  #            overwrite=TRUE)
  #save(elevations,file=paste0(datadir,"/",statesInMap[[1]],"elevs.rda"))
  nchunks <- ceiling(raster::ncell(elevations)/maxrastercells)
  if (nchunks == 1) {
    writeRaster(elevations,file=paste0(datadir,"/rasterfiles/",
                                       statesInMap[[1]],"elevs.grd"),
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
                  file=paste0(datadir,"/rasterfiles/",
                              statesInMap[[1]],"elevs",
                              stringr::str_pad(chunk,2,pad="0"),".grd"),
                  overwrite=TRUE)    
    }
  }
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
                                aspectratio = list(x=1,y=yscale,z=0.03*vertscale),
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
    rgl::aspect3d(x=1,y=1/yscale,z=0.035*vertscale)
    rgl::rgl.clear("lights")
    rgl::rgl.light(theta = 0, phi = 15,
                   viewpoint.rel=TRUE, specular="black")
    rgl::rgl.viewpoint(userMatrix=userMatrix,type="modelviewpoint")
    pan3d(2)  # right button for panning, doesn't play well with zoom)
    if (saveRGL) rgl::writeWebGL(dir=paste0(mapoutputdir), filename=paste0(mapoutputdir,"/",outputName," rgl map.html"))
  }
}
