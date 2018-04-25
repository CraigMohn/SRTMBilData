expandRegions <- function(invec,country) {
  
  ##  sub in vec of states for matching region in apprpriate country
  regionlist <- c("PacificNorthWest",
                  "MountainWest",
                  "DelMarVa",
                  "NewEngland",
                  "Maritimes")
  regionstates <- list(c("WA","OR","ID","MT"),
                       c("WA","OR","ID","MT","WY","CO","UT","NV"),
                       c("DE","MD","DC","VA"),
                       c("ME","NH","VT","MA","RI","CT"),
                       c("NS","PE","NB","NF")
  )
  regioncountry <- c("US","US","US","US","CANADA")
  
  outvec <- invec
  for (i in 1:length(regionlist)) {
    if ((toupper(regionlist[[i]]) %in% toupper(outvec)) &
        (toupper(regioncountry[[i]]) == toupper(country))) {
      outvec <- setdiff(outvec,toupper(regionlist[[i]]))
      outvec <- union(outvec,toupper(regionstates[[i]]))
    }
  }
  return(outvec)
}  

# cropbox <- raster::extent(-160.25, -154.8, -18.9, 22.25) # hawaii main islands only
# cropbox <- raster::extent(-172,-160.25, -18.9, 40) # hawaii NW only


# mapWindow <- c(-159.90,-159.15,21.75,22.35) # Kauai 
# mapWindow <- c(-156.80,-155.97,20.50,21.05) # Maui
# mapWindow <- c(-158.30,-157.60,21.20,21.75) # Oahu
# mapWindow <- c(-123,-121.6,46.75,48)        # Seattle Area 
# mapWindow <- c(-122.2,-121.7,47.4,47.8)     # Samm Area 
# mapWindow <- c(-122.4,-122.1,47.5,47.8)     # Lake Wash Area 
# mapWindow <- c(-123.2,-122.4,48.3,48.8)     # San Juans
# mapWindow <- c(-122.5,-121.9,37.6,38.1)     # East Bay 
# mapWindow <- c(-122.7,-121.8,37.0,38.3)     # SF Bay 
# mapWindow <- c(-122.45,-121.85,37.05,37.65) # Penninsula and South Bay 
# mapWindow <- c(-123.1,-122.42,37.81,38.25)  # Marin County 
# mapWindow <- c(-122.45,-121,40.3,41.6)      # CA Volcanos 
# mapWindow <- c(-81.4,-80.0,36.8,37.6)       # Giles Cty/Blacksburg Area 




