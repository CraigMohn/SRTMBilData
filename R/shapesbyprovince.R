library(raster)
library(rgeos)
library(sp)

datadir <- "c:/bda"                    #  base dir - subs include shapefiles, rasterfile 

#  download data from census canada, unzip into datadir/CanadaShapefiles/
#     creates shapefiles for province roads, cities, water in directory datadir/shapefiles

########################### can skip to load if restarting from saved results

canada <- raster::getData("GADM",country="CAN",level=1) # spatial dataframe

###   cities and towns spdf
tmpt <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/MetroAreas/",
                                 "lcma000b16a_e.shp"))
#  select from CMATYPE B=Metro, K=agglomeration w/tracts, D=agglomeration w/o tracts
tmpt <- tmpt[tmpt@data[,"CMATYPE"]=="B",]

###   highways sldf
tmpr <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/Roads/",
                                 "lrnf000r16a_e.shp"))
#  keep based on RANK 1=transcanada, 2=national, 3=major, 4=secondary hwy
tmpr <- tmpr[as.numeric(tmpr@data[,"RANK"])<=4,]
#  and CLASS - 10-13 hwys, 20-22 roads down to "collector", 23-26 minor roads, 27 rapid transit,
#              80 - bridge/tunnel 
tmpr <- tmpr[as.numeric(tmpr@data[,"CLASS"])<=13,]
rdata <- tmpr@data
row.names(rdata) <- rdata$NGDUID
tmpr <- rgeos::gLineMerge(tmpr,byid=TRUE,id=tmpr@data$NGDUID)
tmpr <- sp::SpatialLinesDataFrame(tmpr,data=rdata)
rdata <- NULL

### areal water spdf
tmpa <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/Lakes+RiversPolygons/",
                                 "lhy_000c16a_e.shp")) 
### linear water sldf
tmpl <- raster::shapefile(paste0(datadir,"/CanadaShapefiles/RiversLines/",
                                 "lhy_000d16a_e.shp")) 
ldata <- tmpl@data
row.names(ldata) <- ldata$HYDROUID
tmpl <- rgeos::gLineMerge(tmpl,byid=TRUE,id=tmpl@data$HYDROUID)
tmpl <- sp::SpatialLinesDataFrame(tmpl,data=ldata)
ldata <- NULL

tmpt <- sp::spTransform(tmpt,crs(canada)) 
tmpr <- sp::spTransform(tmpr,crs(canada)) 
tmpa <- sp::spTransform(tmpa,crs(canada)) 
tmpl <- sp::spTransform(tmpl,crs(canada)) 

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
  tmpmasked <- tmpt[keep,]
  plot(tmpmasked, col="PeachPuff")
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"Town.shp"),
                    overwrite=TRUE)
  print("starting roads")
  keep <- tmpr@data$PRUID_L == prid | tmpr@data$PRUID_L == prid
  tmpmasked <- tmpr[keep,]
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"Roads.shp"),
                    overwrite=TRUE)
  print("starting water area")
  keep <- tmpa@data$PRUID == prid 
  tmpmasked <- tmpa[keep,]
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"WaterA.shp"),
                    overwrite=TRUE)
  print("starting water lines")
  keep <- tmpl@data$PRUID == prid 
  tmpmasked <- tmpl[keep,]
  plot(tmpmasked)
  raster::shapefile(tmpmasked,filename=paste0(datadir,"/shapefiles/",
                                              pr,"WaterL.shp"),
                    overwrite=TRUE)
}


