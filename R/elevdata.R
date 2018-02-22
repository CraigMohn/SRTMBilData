maplist <- c("WA","HI","OR","CA",
                  "ID","UT","NV","AZ","NM",
                  "CO","MN","DelMarVA+NC+WV","WY","MT")
statelist <- list("WA","HI","OR","CA",
                  "ID","UT","NV","AZ","NM",
                  "CO","MN",c("DE","MD","VA","NC","WV"),"WY","MT")
bestmaxelev <- c(3000,2000,3000,3000,
                 3000,3500,3000,3000,3500,
                 3500,1500,1500,3500,3500)

stateregion <- "ID"
res3dplot <- 1300
trimraster <- TRUE

datadir <- "c:/bda"                    #  trimmed elevation raster data location
statedatadir <- "c:/bda/states"        #  zip file subdirectories
mapoutputdir <- "c:/bda/maps3d"        #  map file outputs 
#################################################################################
assign("last.warning", NULL, envir = baseenv())



if (sum(maplist==stateregion)!=1) stop("map state/region error")
highelevation <- bestmaxelev[[which.max(maplist==stateregion)]]
stateMask <- ifelse(trimraster,statelist[[which.max(maplist==stateregion)]],"")
state <- stateregion

# put the srtm data zip files in directory datadir/state
# maps and .rda file with elevations are put in datadir

library(tidyverse)
library(raster)
library(plotly)
library(htmlwidgets)
library(rgl)
library(sp)
library(rgdal)
library(rgeos)
library(tigris)

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


if (stateMask != "") {
  #state <- rgdal::readOGR(dsn = path.data, layer = "usa_state_shapefile")
  #projection(state) <- rgdal::CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
  #state.sub <- state[as.character(state@data$STATE_NAME) %in% stateMask, ]
  state.sub <- rgeos::gUnaryUnion(tigris::counties(stateMask))
  projection(state.sub) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs")
  plot(state.sub)
}

tempd <- tempdir()
r.list <- list()
fn <- list.files(path=paste0(statedatadir,"/",state),pattern="*.zip")
rn <- gsub("_bil","",tools::file_path_sans_ext(fn))
j <- 1
for (i in 1:length(fn)) {
  print(paste0(statedatadir,"/",state,"/",fn[[i]]))
  unzip(paste0(statedatadir,"/",state,"/",fn[[i]]),exdir=tempd)
  tmp <- raster(paste0(tempd,"/",rn[[i]],".bil"))
  if (stateMask != "") {
    xmin <- tmp@extent@xmin
    xmax <- tmp@extent@xmax
    ymin <- tmp@extent@ymin
    ymax <- tmp@extent@ymax
    pgon <- sp::Polygon(cbind(c(xmin,xmax,xmax,xmin,xmin),
                              c(ymin,ymin,ymax,ymax,ymin)))
    ei <- sp::SpatialPolygons(list(Polygons(list(pgon), ID = "1deg")),
             proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
    #ei <- as(extent(tmp), "SpatialPolygons",
    #         proj4string=CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"))
    if (rgeos::gContainsProperly(state.sub, ei)) {
      print("interior - not masked")
    } else if (gIntersects(state.sub, ei)) {
      print(system.time(
               tmp <- raster::mask(raster::crop(tmp, extent(state.sub)),
                                   state.sub)
                        )[[3]])
    } else {
      print("exterior - not used")
      tmp <- NULL
    }
    print(warnings())
  }
  if (!is.null(tmp)) {
    r.list[[j]] <- tmp
    j <- j + 1
  }
}
m.sub <- do.call(merge, r.list)

elevations <- readAll(m.sub)  #  pull it all into memory to avoid complications
save(elevations,file=paste0(datadir,"/",state,"elevs.rda"))
#writeRaster(elevations,file=paste0(datadir,"/",state,"elevs.grd"))

print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- max(1,floor(max(elevations@ncols,elevations@nrows)/res3dplot))
print(paste0("scaling raster down by a factor of ",sfact))
if (sfact > 1)
  elevations <- raster::aggregate(elevations,fact=sfact,fun=mean,
                           expand=TRUE,na.rm=TRUE)
print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))


yscale <- yRatio(elevations)
elevations[is.na(elevations)] <- 0
mmm <- raster::as.matrix(elevations)
mmm <- mmm[,ncol(mmm):1]  #  flip east/west since row 1 is top

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
                            aspectratio = list(x=1,y=yscale,z=0.03),
                            camera=list(up=c(0,1,0),
                                        eye=c(0,1.25,0)) ) )
p
htmlwidgets::saveWidget(p,paste0(mapoutputdir,"/",state," map.html"))

mmmrgl <- raster::as.matrix(elevations)
x <- seq(1,length.out=nrow(mmmrgl))
y <- seq(1,length.out=ncol(mmmrgl))

terrcolors <- colorRampPalette(c("blue","turquoise","aquamarine",
                                 "palegreen","yellowgreen",
                                 "chartreuse","greenyellow","green",
                                 "limegreen","forestgreen","darkgreen",
                                 "yellow","gold","goldenrod",
                                 "sienna","brown","gray75",
                                 "gray85","gray95","gray97","white"))(201)
colidx <- floor(200*(mmmrgl/highelevation)) + 1
colidx[colidx>201] <- 201
colidx[colidx<1] <- 1
col <- terrcolors[colidx]

rgl::par3d("windowRect"= c(100,100,1200,1000))
userMatrix <- matrix(c(-0.02,-0.80,0.632,0,1,0,0.04,0,
                       -0.03,0.60,0.80,0,0,0,0,1),ncol=4,nrow=4)
rgl::rgl.clear()
rgl::surface3d(x,y,mmmrgl,color=col)
rgl::material3d(alpha=1.0,point_antialias=TRUE,smooth=TRUE,shininess=0)
rgl::aspect3d(x=1,y=1/yscale,z=0.02)
rgl::rgl.clear("lights")
rgl::rgl.light(theta = 0, phi = 25,
               viewpoint.rel=TRUE, specular="black")
rgl::rgl.viewpoint(userMatrix=userMatrix,type="modelviewpoint")
pan3d(2)  # right button for panning, doesn't play well with zoom)
rgl::writeWebGL(dir=paste0(mapoutputdir), filename=paste0(mapoutputdir,"/",state," rgl map.html"))

