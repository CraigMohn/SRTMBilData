# SRTMBilData
Script to unzip SRTM .bil files in a map data directory, read and join the raster objects, and save in a large .grd file.  Also processes US TIGER spatial data on urban areas, roads, and water.  The data can be saved in R's native raster .grd format for future use, and can be used from disk rather than being read into memory.  These files are very compressible, and I recommend storing them on a disk where the OS handles that.  Data access is unlikely to be the performance bottleneck - the rgeos functions called in this code can be very slow.

The functions here can also plot a 3d surface and paths on it.  The left mouse button rotates the surface, the right button pans, and the wheel zooms.  Optinally, the function saves an html file with the plot embedded.  The saved plot lacks right mouse-button panning control.  Files saved at a resolution of 2800 (this is a factor used to calculate an aggregation factor for the spatial data) are able to be loaded into Firefox but not Edge or Chrome on my 32G windows machine at this moment, this will likely change with browser releases.

This script is not as graceful as I would like near the 180E/180W meridian.  Fortunately, pretty much nobody lives there.  If you live in the Aleutians, Fiji, Siberia or Chatham Island, I'm sorry for the inconvenience.....

Shuttle Radar Tomography Mission data containing terrestrial surface elevation data is available in several combinations of resolution and void-fill post-processing.  You can create an account and bulk download the files at https://earthexplorer.usgs.gov/  There are R packages which will retrieve SRTM tiles, but the data has been around for a while, and keeping a local copy seems more sensible

There is also a script to create road/city/water shapefiles from Statistics Canada data.

Better instructions to come soon.