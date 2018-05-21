library(raster)
library(rgeos)
library(sp)

#  script to extract Canada shapefiles for towns, roads and water, do rough filtering, 
#      format to match US versions, and save in directory

datadir <- "c:/bda"                    #  base dir - subs include shapefiles, rasterfile 
workProj4 <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs"


###########################################################################################
#  download data from census canada, unzip into datadir/CanadaShapefiles/
#     creates shapefiles for province roads, cities, water in directory datadir/shapefiles
########################### can skip to load if restarting from saved results

canada <- raster::getData("GADM",country="CAN",level=1) # spatial dataframe
canada <- sp::spTransform(canada,workProj4)

###   cities and towns spdf
tmpt <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/MetroAreas/",
                                 "lcma000b16a_e.shp"))
tmpt <- sp::spTransform(tmpt,workProj4)
# type is   B=Metro, K=agglomeration w/tracts, D=agglomeration w/o tracts
tmpt <- tmpt[tmpt@data[,"CMATYPE"]!="K",]

###   highways sldf
tmpr <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/Roads/",
                                 "lrnf000r16a_e.shp"))
#  too huge - keep based on class
#   10 Highway, 11 Expressway, 12 Primary highway, 13 Secondary highway
#   20 Road, 21 Arterial, 22 Collector, 23 Local, 24 Alley/Lane/Utility
#   25 Connector/Ramp, 26 Reserve/Trail, 27 Rapid transit
#   80 - bridge/tunnel 
tmpr <- tmpr[(as.numeric(tmpr@data[,"CLASS"])<=27 | 
              as.numeric(tmpr@data[,"CLASS"])==80) , ]
rdata <- tmpr@data
row.names(rdata) <- rdata$NGDUID
tmpr <- rgeos::gLineMerge(tmpr,byid=TRUE,id=tmpr@data$NGDUID)
tmpr <- sp::SpatialLinesDataFrame(tmpr,data=rdata)

### areal water spdf
tmpa <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/Lakes+RiversPolygons/",
                                 "lhy_000c16a_e.shp")) 
colnames(tmpa@data) <- sub("FULLNAME","NAME",colnames(tmpa@data))

### linear water sldf
tmpl <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/RiversLines/",
                                 "lhy_000d16a_e.shp")) 
ldata <- tmpl@data
row.names(ldata) <- ldata$HYDROUID
colnames(ldata) <- sub("FULLNAME","NAME",colnames(ldata))
tmpl <- rgeos::gLineMerge(tmpl,byid=TRUE,id=tmpl@data$HYDROUID)
tmpl <- sp::SpatialLinesDataFrame(tmpl,data=ldata)

tmpt <- sp::spTransform(tmpt,workProj4) 
tmpr <- sp::spTransform(tmpr,workProj4) 
tmpa <- sp::spTransform(tmpa,workProj4) 
tmpl <- sp::spTransform(tmpl,workProj4) 

save(list = ls(all.names = TRUE), 
     file = paste0(datadir,"/CanadaShapefiles/shapesb4masking.RData"), 
     envir = .GlobalEnv)
# load(file = paste0(datadir,"/CanadaShapefiles/shapesb4masking.RData"))

provlist <- c("BC","AB","SK","MB","ON","QC","NB","NS","NF","PE")
provcode <- c("59","48","47","46","35","24","13","12","10","11")

for (i in 1:length(provlist)) {
  pr <- provlist[[i]]
  prid <- provcode[[i]]
  print(pr)
  mcrop <- rgeos::gUnaryUnion(canada[canada$HASC_1 %in% paste0("CA.",pr),])
  plot(mcrop,col="magenta")
  
  print("starting towns")
  keep <- tmpt@data$PRUID == prid 
  tmpmasked <- tmpt[keep,(names(tmpt) %in% c("CMANAME","CMATYPE"))]
  colnames(tmpmasked@data) <- c("NAME","TYPE")
  plot(tmpmasked, col="PeachPuff")
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"Town.shp"),
                    overwrite=TRUE)
  print("starting roads")
  keep <- tmpr@data$PRUID_L == prid | tmpr@data$PRUID_L == prid
  tmpmasked <- tmpr[keep,(names(tmpr) %in% c("NAME","CLASS"))]
  colnames(tmpmasked@data) <- c("NAME","TYPE")  
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"Roads.shp"),
                    overwrite=TRUE)
  print("starting water area")
  keep <- tmpa@data$PRUID == prid 
  tmpmasked <- tmpa[keep,"NAME"]
  tmpmasked@data$TYPE <- "CANADA"
  tmpsize <- raster::area(tmpmasked)
  tmpsize[tmpsize>10000000] <- 10000000
  tmpmasked@data$size <- tmpsize
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"WaterA.shp"),
                    overwrite=TRUE)
  print("starting water lines")
  keep <- tmpl@data$PRUID == prid 
  tmpmasked <- tmpl[keep,"NAME"]
  tmpmasked@data$TYPE <- "CANADA"
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"WaterL.shp"),
                    overwrite=TRUE)
}


