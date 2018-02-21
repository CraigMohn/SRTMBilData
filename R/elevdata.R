state <- "WA"
datadir <- "c:/bda"
highelevation <- 3000

# put the srtm data zip files in directory datadir/state
# maps and .rda file with elevations are put in datadir

library(tidyverse)
library(raster)
library(plotly)
library(htmlwidgets)
library(rgl)

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


tempd <- tempdir()
r.list <- list()
fn <- list.files(path=paste0(datadir,"/",state),pattern="*.zip")
rn <- gsub("_bil","",tools::file_path_sans_ext(fn))
for (i in 1:length(fn)) {
  print(paste0(datadir,"/",state,"/",fn[[i]]))
  unzip(paste0(datadir,"/",state,"/",fn[[i]]),exdir=tempd)
  r.list[[i]] <- raster(paste0(tempd,"/",rn[[i]],".bil"))
}
m <- do.call(merge, r.list)

elevations <- readAll(m)  #  pull it all into memory to avoid complications
save(elevations,file=paste0(datadir,"/",state,"elevs.rda"))

print(paste0(elevations@ncols," columns by ",elevations@nrows," rows"))
sfact <- max(1,floor(max(elevations@ncols,elevations@nrows)/1500))
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
htmlwidgets::saveWidget(p,paste0(datadir,"/",state," map.html"))

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
rgl::writeWebGL(dir=paste0(datadir), filename=paste0(datadir,"/",state," rgl map.html"))

