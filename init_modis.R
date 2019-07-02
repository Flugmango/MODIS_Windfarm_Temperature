#init_modis

library(devtools)

#install dependencies

# install_github('loicdtx/bfastSpatial')
# install_github("MLezamaValdes/LocST") 
# install.packages("greenbrown", repos="http://R-Forge.R-project.org")

#load dependencies

library(LocST)
library(bfastSpatial)
library(greenbrown)
library(rgdal)
library(MODIS)
library(mapview)
library(bfastSpatial)
# library(Kendall)
# library(ggplot2)
library(latticeExtra) # for as.layer function
library(gridExtra) # for grid arrangement of spplots
# library(rts)
# test <- apply.yearly(rc_full_lst_final, mean)
#preparation
setwd("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature")
lap = "./MODIS"
odp = file.path(lap, "PROCESSED")
MODISoptions(localArcPath = lap, outDirPath = odp)

#adding wind mill positions
windmill_kml <- "Waterloo_Windmills.kml"
windmill_locs <- readOGR(windmill_kml, "Waterloo_Windmills")
slotNames(windmill_locs)
coords <- slot(windmill_locs, "coords")[,-3]
windmill_points <- SpatialPoints(coords)
#----------------------------
#WGS84 APPROACH
# terra orbit tracks https://www.ssec.wisc.edu/datacenter/terra/GLOBAL.html
# bbox creation
bbox= extent(138.858948, 139.113693, -34.355908, -33.966142) #small
bbox2= extent(138.578796, 139.273682, -34.437496, -33.664925) #big
# TERRA LST
# 2003-2010 18 missing   No MOD11A1 files found for the period from 2003-12-17 to 2003-12-23 (351-357):  No MOD11A1 files found for the period from 2008-12-21 to 2008-12-22 (356 - 357)
# 2010-2018 18 missing  No MOD11A1 files found for the period from 2016-02-19 to 2016-02-27 (050 - 058)
#Download Daat
runGdal( product="MYD11A1", extent = bbox2, begin="2000001", end="2018365", SDSstring="10001", outProj="4326") # alternative: begin=as.Date("2000-1-1"), end=as.Date("2018-12-31")

# writeRaster(files_stack, filename="multilayer.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

# #bfast Analysis
# # 5 minutes
# bfast_stack <- timeStackMODIS(files)
# bfast_stack_aqua <- timeStackMODIS(aqua_files)
# # 1.5 hours
# bfast_stack_celsius <- bfast_stack / 50 - 273.15
# bfast_stack_celsius_aqua <- bfast_stack_aqua / 50 - 273.15
# 
# save(bfast_stack_celsius, file="./temp/terra_bfast_stack_celsius.RData")
# load("./temp/terra_bfast_stack_celsius.RData")
# 
# #reclassify files for frost-/no-frost-events
# # reclassify the values into two groups: below 1°C = 1, higher than 1°C = 0
# rc <- reclassify(bfast_stack_celsius, c(-Inf,1,1,1,Inf,0))
# rc_aqua <- reclassify(bfast_stack_celsius_aqua, c(-Inf,1,1,1,Inf,0))
# count_frost <- countObs(rc, navalues = c(0, NA))
# count_frost_aqua <- countObs(rc_aqua, navalues = c(0, NA))
# summary_stack <- summaryBrick(rc, fun=mean, na.rm=TRUE)
# save(summary_stack, file="./temp/summary_stack.RData")
# load("./temp/summary_stack.RData")
# mapview(summary_stack)
# mapview(count_frost)
# 
# # create pre and after aqua data
# aqua_early <- list.files(path="./MODIS/PROCESSED/AQUA_2003_2010", pattern="*.tif$", full.names=TRUE)
# aqua_late <- list.files(path="./MODIS/PROCESSED/AQUA_2010_2018", pattern="*.tif$", full.names=TRUE)
# aqua_early_stack <- timeStackMODIS(aqua_early)
# aqua_late_stack <- timeStackMODIS(aqua_late)
# #change reclassify
# # rc_aqua_early <- reclassify(aqua_early_stack, c(-Inf,1,1,1,Inf,0))
# rc_aqua_early <- reclassify(aqua_early_stack, c(-Inf,13707.5,1,13707.5,Inf,0))
# save(rc_aqua_early, file="./temp/rc_aqua_early.RData")
# load("./temp/rc_aqua_early.RData")
# rc_aqua_late <- reclassify(aqua_late_stack, c(-Inf,13707.5,1,13707.5,Inf,0))
# save(rc_aqua_late, file="./temp/rc_aqua_late.RData")
# load("./temp/rc_aqua_late.RData")
# frost_aqua_early <- countObs(rc_aqua_early, navalues = c(0, NA))
# frost_aqua_late <- countObs(rc_aqua_late, navalues = c(0, NA))
# aqua_frost_stack <- stack(c(frost_aqua_early, frost_aqua_late))
# aqua_names <- c("aqua_2003_2010", "aqua_2010_2018")
# names(aqua_frost_stack) <- aqua_names
# spplot(aqua_frost_stack)
# 
# # create pre and after terra data
# terra_early <- list.files(path="./MODIS/PROCESSED/TERRA_2003_2010", pattern="*.tif$", full.names=TRUE)
# terra_late <- list.files(path="./MODIS/PROCESSED/TERRA_2010_2018", pattern="*.tif$", full.names=TRUE)
# terra_early_stack <- timeStackMODIS(terra_early)
# terra_late_stack <- timeStackMODIS(terra_late)
# 
# #change reclassify
# # rc_aqua_early <- reclassify(aqua_early_stack, c(-Inf,1,1,1,Inf,0))
# rc_terra_early <- reclassify(terra_early_stack, c(-Inf,13707.5,1,13707.5,Inf,0))
# save(rc_terra_early, file="./temp/rc_terra_early.RData")
# load("./temp/rc_terra_early.RData")
# rc_terra_late <- reclassify(terra_late_stack, c(-Inf,13707.5,1,13707.5,Inf,0))
# save(rc_terra_late, file="./temp/rc_terra_late.RData")
# load("./temp/rc_terra_late.RData")
# full_rc_terra <- stack(rc_terra_early, rc_terra_late)
# frost_terra_early <- countObs(rc_terra_early, navalues = c(0, NA))
# frost_terra_late <- countObs(rc_terra_late, navalues = c(0, NA))
# terra_frost_stack <- stack(c(frost_terra_early, frost_terra_late))
# terra_names <- c("terra_2003_2010", "terra_2010_2018")
# names(terra_frost_stack) <- terra_names
# spplot(terra_frost_stack)
# 
# ### NUMBER OF FROST EVENT PER YEAR ACCUMULATION
# 
# 
# 
# #get the date from the names of the layers and extract the year
# terra_2003_2010 <- list.files(path="./MODIS/PROCESSED/TERRA_2003_2010", pattern="*.tif$", full.names=TRUE)
# foo <- timeStackMODIS(terra_2003_2010)
# years <- getZ(foo)
# names(foo) <- years
# #indices <- format(as.Date(names(full_terra), format = "X%Y.%m.%d"), format = "%Y")
# indices <- format(as.Date(getZ(foo), format = "X%Y.%m.%d"), format = "%Y")
# indices <- as.numeric(indices)
# # ~2 hours for full_terra
# rc_full_terra <- reclassify(full_terra, c(-Inf,13707.5,1,13707.5,Inf,0))
# save(rc_full_terra, file="./temp/rc_full_terra.RData")
# load("./temp/rc_full_terra.RData")
# #maybe setZ(rc_full_terra, years, "time")
# setZ(rc_full_terra, as.Date(years), "time")
# getZ(rc_full_terra)

################ MAIN
#Loading data
terra_lst_files <- list.files(path="./MODIS/PROCESSED/MODIS_TERRA_LST_2000_2018", pattern="*.tif$", full.names=TRUE)
aqua_lst_files <- list.files(path="./MODIS/PROCESSED/MODIS_AQUA_LST_2002_2018", pattern="*.tif$", full.names=TRUE)
# processing data
# ~ 8mins
start_time <- Sys.time()
stack = timeStackMODIS(c(terra_lst_files, aqua_lst_files))
end_time <- Sys.time()
save(stack, file="./temp/stack.RData")
load("./temp/stack.RData")
#order stack in tiem
new_order <- order(getZ(stack))
stack <- subset(stack, new_order)
# names(full_terra) <- years
years <- getZ(stack)
# ~5 hours for full
rc_full_lst <- reclassify(stack, c(0,13707.5,1,13707.5,65535,0))
save(rc_full_lst, file="./temp/rc_full_lst.RData")
load("./temp/rc_full_lst.RData")
#---------------------------
# ~5 hours
rc_full_lst_ordered <- reclassify(stack, c(0,13707.5,1,13707.5,65535,0))
save(rc_full_lst_ordered, file="./temp/rc_full_lst_ordered.RData")
load("./temp/rc_full_lst_ordered.RData")
rc_full_lst_final <- setZ(rc_full_lst_ordered, years)

start_time <- Sys.time()
trendmap <- TrendRaster(rc_full_lst_final, start=c(2000, 2), freq=365, method="AAT", breaks=2)
end_time <- Sys.time()
spplot(trendmap2) + as.layer(spplot(windmill_points, pch = 3, edge.col = "white"))


trendmap2 <- TrendRaster(rc_full_lst_final)
plot(trendmap2)

#experimental
plots <- lapply(names(trendmap2[1:3]), function(.x) spplot(trendmap2, .x))
do.call(grid.arrange, plots)
grid.arrange(spplot(df.voro2, "my.data", xlab = "x", ylab = "y", main = "my title") , 
             spplot(df.voro, "dummy", xlab = "x2", ylab = "y2", main = "my title2" ))
#---------------------------
count_full_lst <- countObs(rc_full_lst_final, navalues = c(0, NA))
spplot(count_full_lst) + as.layer(spplot(windmill_points, pch = 3, edge.col = "white"))

#build mean of the layers
indices <- format(as.Date(years, format = "X%Y.%m.%d"), format = "%Y")
indices <- as.numeric(indices)
datasum<- stackApply(rc_full_lst_final, indices, fun = mean)
spplot(datasum)

#get prettier names for plot titles
datasum_names <- names(datasum)
datasum_names <- strsplit(datasum_names, "_")[[1]][2]

for (i in 1:length(datasum_names)) {
    datasum_names[i] <- strsplit(datasum_names[i], "_")[[1]][2]
}

names(datasum) <- datasum_names
substr(names(datasum), 2, 5)
spplot(datasum)

####### extract exact time stamp

# ltz <- "Australia/Adelaide"   # https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
# newproj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

lst <- raster("./MODIS/PROCESSED/MOD11A1_VT_2001_2018/MOD11A1.A2001001.Day_view_time.tif")
utc <- strptime(c("2001001"), format="%Y%j", tz="UTC") #get utc day from filename with strptime(c("2016001"), format="%Y%j", tz="UTC")
ltz <- "Australia/Adelaide"

locst_results <- LocST_UTC(lst, utc)
utc_times <- locst_results$utc_from_LocST
acdt_times <- format(as.POSIXct(utc_times), usetz= TRUE, tz="Australia/Adelaide")

zApply()
# using lmvt
fnam <- "./MODIS/PROCESSED/MOD11A1.006_20190605184217/MOD11A1.A2017055.LST_Day_1km.tif"
viewtime <- "./MODIS/PROCESSED/MOD11A1.006_20190605184217/MOD11A1.A2017055.Day_view_time.tif"
x <- as.POSIXct(locst_results$utc_from_LocST[1])
y <- format(x, tz="Australia/Adelaide", usetz=TRUE)

lmvt(fnam = lst, viewtime = utc, ltz = ltz)






#Kendall analysis

# fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
fun_kendall <- function(x){ return(unlist(MannKendall(x)))}
kendall_result <- calc(celsius_stack,fun_kendall)

x1 <- writeRaster(kendall_result, filename = "./kendall_output/MODIS_terra_layer_1.tif", format="GTiff", overwrite=TRUE)
x2 <- writeRaster(kendall_result$layer.2, filename = "./kendall_output/MODIS_terra_layer_2.tif", format="GTiff", overwrite=TRUE)
x3 <- writeRaster(kendall_result$layer.3, filename = "./kendall_output/MODIS_terra_layer_3.tif", format="GTiff", overwrite=TRUE)
x4 <- writeRaster(kendall_result$layer.4, filename = "./kendall_output/MODIS_terra_layer_4.tif", format="GTiff", overwrite=TRUE)
x5 <- writeRaster(kendall_result$layer.5, filename = "./kendall_output/MODIS_terra_layer_5.tif", format="GTiff", overwrite=TRUE)
