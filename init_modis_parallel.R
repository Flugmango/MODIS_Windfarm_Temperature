#init_modis_parallel

#load dependencies
library(rgdal)
library(MODIS)
library(mapview)
library(parallel)

#prepare parallelization
#calculate number of cores
no_cores = max(1, detectCores() -1)
#initiate cluster
cl = makeCluster(no_cores)

setwd("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature")
#preparation
lap = "./MODIS"
odp = file.path(lap, "PROCESSED")
MODISoptions(localArcPath = lap, outDirPath = odp)

#----------------------------
#WGS84 APPROACH
# bbox creation
bbox= extent(138.858948, 139.113693, -34.355908, -33.966142)
#Download Data
runGdal( product="MOD11A1", extent = bbox, begin="2018001", end="2018365", SDSstring="1010101", outProj="4326")

# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000001", end="2018365", SDSstring="10001", outProj="4326")
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin=as.Date("2000-1-1"), end=as.Date("2018-12-31"), SDSstring="10001", outProj="4326")

# processing data
files <- list.files(path="./MODIS/PROCESSED/MOD11A1.006_20190314172848", pattern="*.tif$", full.names=TRUE)
stack = stack(files)

# save(stack, file="./temp/terra_stack.RData")
# load("./temp/terra_stack.RData")

# clusterExport(cl, "files") #not sure if needed
parLapply(cl, stack, function(raster_kelvin){
  raster_celsius = (raster_kelvin / 50) - 273.15
  # save(raster_celsius, file="files_celsius.RData")
}
)
# writeRaster(files_stack, filename="multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
# downscale and transform kelvin
#takes almost 45 minutes
celsius_stack = (stack / 50) - 273.15
#saving
save(celsius_stack, file="./temp/terra_stack_celsius.RData")
load("./temp/terra_stack_celsius.RData")
#free the resources
stopCluster(cl)








# function for processing rasters
# parLapply(files, function(x) {
#   #read .tif
#   curr_raster <- readGDAL(x)
#   # crop
#   out <- function(curr_raster) {
#   }
#     # write to file
#     write.table(out, "./MODIS/", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# })

# x <- readGDAL("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature/MODIS/PROCESSED/MOD11A1.006_20190312151749/MOD11A1.A2000056.LST_Day_1km.tif")
# make x croppable
# x = raster(x)
# result = crop(x, bbox)

# mapview(result)
