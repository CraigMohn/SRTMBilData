state <- "WA"
datadir <- "c:/bda"

library(tidyverse)
library(raster)
library(plotly)
library(htmlwidgets)

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
elevations[is.na(elevations)] <- -100
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
                            aspectratio = list(x=1,y=1/1.4,z=0.03),
                            camera=list(up=c(0,1,0),
                                        eye=c(0,1.25,0)) ) )
p
htmlwidgets::saveWidget(p,paste0(datadir,"/",state," map.html"))


