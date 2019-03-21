#init_modis

#load dependencies
library(rgdal)
library(MODIS)
library(raster)
library(mapview)
# library(parallel)

# setwd("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature")
#preparation
lap = "./MODIS"
odp = file.path(lap, "PROCESSED")
MODISoptions(localArcPath = lap, outDirPath = odp)


#----------------------------
#WGS84 APPROACH
# bbox creation
bbox= extent(138.858948, 139.113693, -34.355908, -33.966142)
#Download Data
runGdal( product="MOD11A1", extent = bbox, begin="2011068", end="2018365", SDSstring="10001", outProj="4326")
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000001", end="2018365", SDSstring="10001", outProj="4326")
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin=as.Date("2000-1-1"), end=as.Date("2018-12-31"), SDSstring="10001", outProj="4326")
# x <- readGDAL("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature/MODIS/PROCESSED/MOD11A1.006_20190312151749/MOD11A1.A2000056.LST_Day_1km.tif")
#make x croppable
# x = raster(x)

# processing data
files <- list.files(path="./MODIS/PROCESSED/MOD11A1.006_20190314172848", pattern="*.tif$", full.names=TRUE)
files_stack <- stack(files)
files_50 = files_stack / 50
avg_files_stack = cellStats(files_stack, mean)
lapply(files, function(x) {
  #read .tif
  curr_raster <- readGDAL(x)
  # crop
  out <- function(curr_raster) {
  }
    # write to file
    write.table(out, "./MODIS/", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
})
# result = crop(x, bbox)

# mapview(result)