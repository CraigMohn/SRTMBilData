# SRTMBilData
Script to unzip SRTM .bil files in a "state" directory, read and join the raster objects, and save in a large .rda file

It also calls plotly to plot and try to save the data in an html file.  It often works. The file will be large, and display clumsily.  I will probably soon paste few lines of code to render the maps with rgl, which can be saved as truly unbelievably large files, but the 3d image displayed is much more stable and easy to manipulate .  

Shuttle Radar Tomography Mission data containing terrestrial surface elevation data is available in several combinations of resolution and void-fill post-processing.  You can create an account and bulk download the files at https://earthexplorer.usgs.gov/

