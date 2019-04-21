#init_modis

#load dependencies
library(rgdal)
library(MODIS)
library(mapview)
library(Kendall)

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

# saving
save(stack, file="./temp/terra_stack.RData")
load("./temp/terra_stack.RData")

# writeRaster(files_stack, filename="multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
# 
# downscale and transform kelvin to celsius: takes almost 45 minutes
celsius_stack = (stack / 50) - 273.15

#saving
save(celsius_stack, file="./temp/terra_stack_celsius.RData")
load("./temp/terra_stack_celsius.RData")

#Kendall analysis

# fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
kendall_result <- calc(celsius_stack,fun_kendall)

x1 <- writeRaster(kendall_result, filename = "./kendall_output/MODIS_terra_layer_1.tif", format="GTiff", overwrite=TRUE)
x2 <- writeRaster(kendall_result$layer.2, filename = "./kendall_output/MODIS_terra_layer_2.tif", format="GTiff", overwrite=TRUE)
x3 <- writeRaster(kendall_result$layer.3, filename = "./kendall_output/MODIS_terra_layer_3.tif", format="GTiff", overwrite=TRUE)
x4 <- writeRaster(kendall_result$layer.4, filename = "./kendall_output/MODIS_terra_layer_4.tif", format="GTiff", overwrite=TRUE)
x5 <- writeRaster(kendall_result$layer.5, filename = "./kendall_output/MODIS_terra_layer_5.tif", format="GTiff", overwrite=TRUE)
