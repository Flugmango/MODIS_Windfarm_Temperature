#init_modis

library(rgdal)
library(MODIS)
library(raster)
library(mapview)

lap = "./MODIS"
odp = file.path(lap, "PROCESSED")
MODISoptions(localArcPath = lap, outDirPath = odp)
# outProj: 4326 = WGS84; 3112 = GDA94 / Geoscience Australia Lambert; 8059 South Australia
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000055", end="2000060", SDSstring="10001", outProj="4326")
runGdal( product="MOD11A1", extent = 'Australia', begin="2000055", end="2000060", SDSstring="10001", outProj="1150")

x <- readGDAL("C:/Users/Boris/Documents/MODIS/PROCESSED/MOD11A1.006_20190217124837/MOD11A1.A2000055.LST_Day_1km.tif")