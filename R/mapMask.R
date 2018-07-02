mapMask <- function(USStatevec=NULL,CAProvincevec=NULL,
                    USParkvec=NULL,worldCountryvec=NULL,
                    mapWindow=NULL,
                    mapbuffer=0,mapmergebuffer=10,
                    parkDir,
                    workProj4="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                    year=2017) {
##    return a spatial polygon mask for the map 
  if (is.null(USStatevec)&is.null(CAProvincevec)&
      is.null(USParkvec)&is.null(worldCountryvec)&is.null(mapWindow))
    stop(paste0("nothing specified for map"))

  mapshape <- NULL
  if (is.null(mapWindow)) {
    if (!is.null(USStatevec) & is.null(USParkvec)) {
      USStatevec <- expandRegions(unique(toupper(USStatevec)),"US")  # US State abbrev all upper case  
      mcrop <- tigris::counties(USStatevec,year=year) %>% 
               rgeos::gUnaryUnion(.) %>% 
               sp::spTransform(.,sp::CRS(workProj4))
      ## tigris returns a SpatialPolgonsDF and gUnaryUnion returns a SpatialPolygons
      mapshape <- bufferUnion(spObj=mcrop,mapbuffer=mapmergebuffer,
                              mapunion=mapshape)
    }
    if (!is.null(USParkvec)) {
      USParkvec <- unique(toupper(USParkvec))
      pfile <- "nps_boundary.shp"
      parkareas <- sf::st_read(paste0(parkDir,"/",pfile))  # sf dataframe
      #  parknames <- parkareas[,c("UNIT_NAME","UNIT_CODE")]
      #  sf::st_geometry(parknames) <- NULL
      #  parknames <- parknames[order(parknames$UNIT_CODE),] 
      #  write.csv(parknames,paste0(datadir,"/parknames.csv"))
      parkareas <-  sf::st_transform(parkareas[parkareas$UNIT_CODE %in% USParkvec,"geometry"],
                                   workProj4) 
      mcrop <- as(parkareas, "Spatial")  
      #  don't include the USStatevec states
      mapshape <- bufferUnion(spObj=mcrop,mapbuffer=mapmergebuffer,
                              mapunion=NULL) 
    }
    if (!is.null(CAProvincevec)) {
      CAProvincevec <- expandRegions(unique(toupper(CAProvincevec)),"CANADA")
      canada <- raster::getData("GADM",country="CAN",level=1) %>%
                sp::spTransform(.,sp::CRS(workProj4))
      mcrop <- canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),] %>%
      rgeos::gUnaryUnion(.)
      # simplify - BC coast is extremely complex
      mapshape <- bufferUnion(spObj=mcrop,mapbuffer=mapmergebuffer,
                              mapunion=mapshape,simplifytol = 1) 
    }
    if (!is.null(worldCountryvec)) { 
      worldCountryvec <- unique(toupper(worldCountryvec))
      mcrop <- NULL
      for (c in worldCountryvec) {
print(c)
        cmap <- raster::getData("GADM",country=c,level=0)
plot(cmap)
        cmap <- rgeos::gUnaryUnion(cmap) 
        cmap <-  sp::spTransform(cmap,sp::CRS(workProj4))
        if (is.null(mcrop)) {
          mcrop <- cmap
        } else {
          mcrop<- rgeos::gUnaryUnion(raster::union(mcrop,cmap))
        }
      }
      mapshape <- bufferUnion(spObj=mcrop,mapbuffer=mapmergebuffer,
                              mapunion=mapshape,simplifytol = 1)
    }
  } else {
    mapshape <- raster::extent(mapWindow)
    CP <- as(mapshape, "SpatialPolygons")
    sp::proj4string(CP) <- workProj4
    mapshape <- rgeos::gUnaryUnion(CP)
  }
  return(bufferUnion(mapshape,mapbuffer=mapbuffer,mapunion=NULL,simplifytol=0))
}