#init_modis
#load dependencies
library(devtools)
# install_github('loicdtx/bfastSpatial')
# install.packages("greenbrown", repos="http://R-Forge.R-project.org")
library(bfastSpatial)
library(greenbrown)
library(rgdal)
library(MODIS)
library(mapview)
library(bfastSpatial)
library(Kendall)

setwd("C:/Users/Flugmango/Documents/MODIS_Windfarm_Temperature")
source("locmodisviewtime_LocST.R")
#preparation
lap = "./MODIS"
odp = file.path(lap, "PROCESSED")
MODISoptions(localArcPath = lap, outDirPath = odp)

#----------------------------
#WGS84 APPROACH
# terra orbit tracks https://www.ssec.wisc.edu/datacenter/terra/GLOBAL.html
# bbox creation
bbox= extent(138.858948, 139.113693, -34.355908, -33.966142)
#Download Data
runGdal( product="MOD11A1", extent = bbox, begin="2001001", end="2018365", SDSstring="1010101", outProj="4326")

# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000001", end="2018365", SDSstring="10001", outProj="4326")
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin=as.Date("2000-1-1"), end=as.Date("2018-12-31"), SDSstring="10001", outProj="4326")

# processing data
files <- list.files(path="./MODIS/PROCESSED/MOD11A1.006_2000_2018", pattern="*.tif$", full.names=TRUE)
# ~ 5mins
stack = stack(files)
end_time <- Sys.time()
#reclassify files for frost-/no-frost-events
# reclassify the values into two groups: below 1°C = 1, higher than 1°C = 0
rc <- reclassify(stack, c(-Inf,1,1, 1,Inf,0))


# saving
save(stack, file="./temp/terra_stack.RData")
load("./temp/terra_stack.RData")

# writeRaster(files_stack, filename="multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
# 
# downscale and transform kelvin to celsius: takes approx. 45 minutes
celsius_stack = (stack / 50) - 273.15
#saving
save(celsius_stack, file="./temp/terra_stack_celsius.RData")
load("./temp/terra_stack_celsius.RData")


#bfast Analysis
# 5 minutes
bfast_stack <- timeStackMODIS(files)
# 1.5 hours
bfast_stack_celsius <- bfast_stack / 50 - 273.15

save(bfast_stack_celsius, file="./temp/terra_bfast_stack_celsius.RData")
load("./temp/terra_bfast_stack_celsius.RData")

summary_stack <- summaryBrick(bfast_stack_celsius, fun=mean, na.rm=TRUE)
save(summary_stack, file="./temp/summary_stack.RData")
load("./temp/summary_stack.RData")
mapview(summary_stack)

#######adding time vector on hold

#converting LST to LT
# fnam <- "./MODIS/MODIS/MOD11A1.006/MOD11A1.A2019137.h29v12.006.2019138085521.hdf"
# ra <- c(138.858948, 139.113693, -34.355908, -33.966142)
# viewtime <- raster("./MODIS/PROCESSED/MOD11A1.006_20190521001332/MOD11A1.A2019137.Day_view_time.tif")
# ltz <- "Australia/Adelaide"   # https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
# newproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# 
# locmodisviewtime(fnam, viewtime, ra, ltz, newproj)


#Kendall analysis

# fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
kendall_result <- calc(celsius_stack,fun_kendall)

x1 <- writeRaster(kendall_result, filename = "./kendall_output/MODIS_terra_layer_1.tif", format="GTiff", overwrite=TRUE)
x2 <- writeRaster(kendall_result$layer.2, filename = "./kendall_output/MODIS_terra_layer_2.tif", format="GTiff", overwrite=TRUE)
x3 <- writeRaster(kendall_result$layer.3, filename = "./kendall_output/MODIS_terra_layer_3.tif", format="GTiff", overwrite=TRUE)
x4 <- writeRaster(kendall_result$layer.4, filename = "./kendall_output/MODIS_terra_layer_4.tif", format="GTiff", overwrite=TRUE)
x5 <- writeRaster(kendall_result$layer.5, filename = "./kendall_output/MODIS_terra_layer_5.tif", format="GTiff", overwrite=TRUE)
