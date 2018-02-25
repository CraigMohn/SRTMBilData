# SRTMBilData
Script to unzip SRTM .bil files in a "state" directory, read and join the raster objects, and save in a large .rda file

It also calls plotly to plot and try to save the data in an html file.  It often works. The file will be large, and display clumsily in the browsers I have tried.  The script also creates an rgl version of the 3d object, and saves an html file with it embedded.  THese files are an order of magnitude larger than the plotly version, and if you can find a browser which successfully opens them, they don't seem to have the spontaneous resets and lockups that the plotly versions exhibit.  Firefox seems to do better with large objects than Edge or Chrome on my machine at this moment.

This script is not as graceful as I would like near the 180E/180W meridian.  Fortunately, nobody lives there.  If you live in the Aleutians, Fiji, Siberia or Chatham Island, I'm sorry for the inconvenience.....

Shuttle Radar Tomography Mission data containing terrestrial surface elevation data is available in several combinations of resolution and void-fill post-processing.  You can create an account and bulk download the files at https://earthexplorer.usgs.gov/

