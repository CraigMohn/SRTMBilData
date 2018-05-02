mapMask <- function(USStatevec=NULL,CAProvincevec=NULL,
                    USParkvec=NULL,worldCountryvec=NULL,
                    mapWindow=NULL,
                    mapbuffer=0,mapmergebuffer=10,
                    parkdir,
                    workProj4="+proj=longlat +ellps=WGS84 +towgs84=0,0,0 +no_defs") {
##    return a spatial polygon mask for the map 
  if (is.null(USStatevec)&is.null(CAProvincevec)&
      is.null(USParkvec)&is.null(worldCountryvec)&is.null(mapWindow))
    stop(paste0("nothing specified for map"))

  mapcrop <- NULL
  if (!is.null(USStatevec)) {
    USStatevec <- expandRegions(unique(toupper(USStatevec)),"US")  # US State abbrev all upper case
    mcrop <- tigris::counties(USStatevec) %>% 
             rgeos::gUnaryUnion(.) %>% 
             sp::spTransform(.,sp::CRS(workProj4))
    ## tigris returns a SpatialPolgonsDF and gUnaryUnion returns a SpatialPolygons
    mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop)
  }
  if (!is.null(CAProvincevec)) {
    CAProvincevec <- expandRegions(unique(toupper(CAProvincevec)),"CANADA")
    canada <- raster::getData("GADM",country="CAN",level=1) %>%
              sp::spTransform(.,sp::CRS(workProj4))
    mcrop <- canada[canada$HASC_1 %in% paste0("CA.",CAProvincevec),] %>%
    rgeos::gUnaryUnion(.)
    # simplify - BC coast is extremely complex
    mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1) 
  }
  if (!is.null(USParkvec)) {
    USParkvec <- unique(toupper(USParkvec))
    pfile <- "nps_boundary.shp"
    parkareas <- sf::st_read(paste0(parkdir,"/",pfile))  # sf dataframe
    #parkareas <- rgeos::gUnaryUnion(parkareas)
    #  parknames <- parkareas[,c("UNIT_NAME","UNIT_CODE")]
    #  sf::st_geometry(parknames) <- NULL
    #  parknames <- parknames[order(parknames$UNIT_CODE),] 
    #  write.csv(parknames,paste0(datadir,"/parknames.csv"))
    parkareas <-  sp::spTransform(parkareas[parkareas$UNIT_CODE %in% USParkvec,"geometry"],
                                  sp::CRS(workProj4)) 
    mcrop <- as(parkareas, "Spatial")  
    mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop) 
  }
  if (!is.null(worldCountryvec)) { 
    worldCountryvec <- unique(toupper(worldCountryvec))
    mcrop <- NULL
    for (c in worldCountryvec) {
      cmap <- raster::getData("GADM",country=worldCountryvec,level=0) # raster + spatial
      cmap <- rgeos::gUnaryUnion(cmap) 
      cmap <-  sp::spTransform(cmap,sp::CRS(workProj4))
      if (is.null(mcrop)) {
        mcrop <- cmap
      } else {
        mcrop<- rgeos::gUnaryUnion(raster::union(mcrop,cmap))
      }
    }
    mapcrop <- bufferUnion(mcrop,mapbuffer=mapmergebuffer,mapcrop,simplifytol = 1)
  }
  ## mapWindow - overwrite the map used for cropping, 
  ##     but not the list of states/provinces to load
  if (!is.null(mapWindow)) {
    mapcrop <- raster::extent(mapWindow)
    CP <- as(mapcrop, "SpatialPolygons")
    sp::proj4string(CP) <- workProj4
    mapcrop <- rgeos::gUnaryUnion(CP)
  }
  return(bufferUnion(mapcrop,mapbuffer=mapbuffer,NULL,simplifytol = 0))
}