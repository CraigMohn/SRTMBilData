filterFeatures <- function(spList,filterList=NULL) {
  #  return a list which is a copy of the input list with each named element
  #    containing the original spatial dataframe with elements less than 
  #    the numeric value with the same name in filterList
  if (is.null(filterList)) 
    filterList <- list("spTown"=1,
                       "spRoads"=1,
                       "spWaterA"=1,
                       "spWaterL"=1)
  
  spListFiltered <- spList
  for (x in names(spList)) {
    if (length(filterList[[x]])>0) {
      temp <- spList[[x]]
      xx <- sub("sp","",x)
      xx <- paste(tolower(substr(xx,1,1)),substr(xx,2,nchar(xx)),sep="")
      rankf <- get(paste0(xx,"Rank"))
      if ("size" %in% names(temp@data)) {
        temp$value <- rankf(temp@data[,"TYPE"],temp@data[,"NAME"],temp@data[,"size"])
      } else {
        temp$value <- rankf(temp@data[,"TYPE"],temp@data[,"NAME"])
      }
    }
    spListFiltered[[x]] <- temp[temp$value>=filterList[[x]],]
  }
  return(spListFiltered)
}
townRank <- function(ttype,tname) {
  ttype <- as.vector(ttype)
  # for US TIGER data type is   U pop >= 50k, 2.5k <= C pop < 50k 
  # for Canada type is   B=Metro, K=agglomeration w/tracts, D=agglomeration w/o tracts
  #    (K is dropped already)
  rankT <- rep(1,length(ttype))
  rankT[ttype %in% c("C","D")] <- 3           
  rankT[ttype %in% c("U","B")] <- 5           
  return(as.integer(rankT))
}
roadsRank <- function(rtype,rname) {
  rtype <- as.vector(rtype)
  # for US TIGER data type 
  # type is   S1100=Primary        S1200=Secondary      S1400=Local St 
  #           S1500=Vehicular Trl  S1630=Ramp           S1640=Service Drive 
  #           S1710=Walkway        S1720=Stairway       S1730=Alley
  #           S1740=Private Rd     S1750=Internal       S1780=Pkgng Lot 
  #           S1820=Bike Path      S1830=Bridle Trl
  # for Canada type is from(CLASS)
  #   10 Highway, 11 Expressway, 12 Primary highway, 13 Secondary highway
  #   20 Road, 21 Arterial, 22 Collector, 23 Local, 24 Alley/Lane/Utility
  #   25 Connector/Ramp, 26 Reserve/Trail, 27 Rapid transit
  #   80 - bridge/tunnel 
  rankR <- rep(1,length(rtype))                    # anything here gets a 1
  rankR[rtype %in% c("S1500","S1630","S1640",
                     "S1730","S1780","S1820",
                     "24","25","26")]         <- 2 # Service Drive, Bike Path, etc
  rankR[rtype %in% c("S1400",
                     "20","21","22","23")]    <- 3 # Local Street
  rankR[rtype %in% c("S1100","10","13","80")] <- 4 # secondary hwy
  rankR[rtype %in% c("S1200","11","12","27")] <- 5 # primary hwy/transit
  return(as.integer(rankR))
}
waterARank <- function(wtype,wname,wsize) {
  wtype <- as.vector(wtype)
  wname <- as.vector(wname)
  wsize <- as.vector(wsize)
  # type is   H2025=Swamp,H2030=Lake/Pond,H2040=Reservoir,H2041=TreatmentPond,
  #           H2051=Bay/Est/Sound,H2053=Ocean,H2060=Pit/Quarry,H2081=Glacier
  # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
  #           CANADA=from canada files
  rankA <- rep(1,length(wtype))                      # anything here gets at least 1
  rankA[wtype %in% c("H2040","H2041","H2060")] <- 2  # res/treatmentpond/pit/quarry
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA")]                <- 3  # lake/pond/swamp/stream
  
  rankA[wtype %in% c("H2040","H2041","H2060",
                     "H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          wsize >= 1000]                        <- 4 # any of the previous that are not small
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          !is.na(wname)]                        <- 5 # lake/pond/swamp/stream named
  rankA[wtype %in% c("H2025","H2030",
                     "H3010","H3013","H3020",
                     "CANADA") &
          wsize >= 4000]                        <- 6 # lake/pond/swamp/stream big
  rankA[wtype %in% c("H2053","H2051")]          <- 7 # Ocean/Bay/Est/Sound
  rankA[wtype == "H2081"]                       <- 8 # glacier
  return(as.integer(rankA))
}
waterLRank <- function(wtype,wname) {
  wtype <- as.vector(wtype)
  wname <- as.vector(wname)
  #  classify the linear water by ranking and sort it so higher rank overwrites lower
  # type is   H3010=Stream/River,H3013=BraidedStream,H3020=Canal/Ditch
  #   H1100 is unspecified, gets  set to 1
  rankL <- rep(1,length(wtype))               # anything here gets a 1
  rankL[wtype == "H3020"] <- 2                # canal/ditch
  rankL[wtype == "H3013"] <- 3                # H3013 braided stream
  rankL[wtype %in% c("H3010","CANADA")] <- 4  # H3010 stream/river
  rankL[wtype %in% c("H3010","CANADA") & 
          !is.na(wname)] <- 5                 # H3010 + name not missing
  rankL[wtype %in% c("H3010","CANADA") &
          grepl("RIV",toupper(wname))] <- 6   # H3010 + name contains "RIV"
  return(as.integer(rankL))
}
