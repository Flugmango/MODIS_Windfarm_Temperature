#backup script

# <-- alternative bbox creation
# bbox_matrix = cbind(c(138.858948, 138.858948, 139.113693, 139.113693, 138.858948),c(-33.966142, -34.355908, -34.355908, -33.966142, -33.966142))
# p = Polygon(bbox_matrix)
# ps = Polygons(list(p),1)
# sps = SpatialPolygons(list(ps))
# proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# temp = data.frame(f=99.9)
# spdf = SpatialPolygonsDataFrame(sps,temp)
# # x = raster(x)
# result = crop(x, spdf)
# -->

#backup
# outProj: 4326 = WGS84; 3112 = GDA94 / Geoscience Australia Lambert; 8059 South Australia
# bbox frost damage GDA94: 443019.7219,-3873566.1926,467169.3498,-3835464.7163
# bbox frost damage WGS84: 138.858948,-34.355908,139.113693,-33.966142
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000055", end="2000060", SDSstring="10001", outProj="28354")
# 28354 +proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs
#MODIS terra data download
# runGdal(product="MOD11A1", extent = 'Australia', begin="2000055", end="2000060", SDSstring="10001", outProj="28354")
# x1 = mapview(bbox_gda94, native.crs= TRUE)
# x2 = mapview(x, native.crs= TRUE)
# latticeView(x1, x2)
# plot(x)

#----------------------------
#CUSTOM CRS APPROACH
# runGdal( product="MOD11A1", tileH = 29, tileV = 12, begin="2000055", end="2000060", SDSstring="10001", outProj="28354")

#Data subsetting
#creating bounding box spatial polygon
# bbox_wgs84 = as(raster::extent(138.858948, 139.113693, -34.355908, -33.966142), "SpatialPolygons")
# proj4string(bbox_wgs84) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# plot(bbox_wgs84)

# bbox_gda94 = spTransform(bbox_wgs84, CRS("+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
# plot(bbox_gda94)

# x <- readGDAL("C:/Users/Boris/Documents/MODIS/PROCESSED/MOD11A1.006_20190217124837/MOD11A1.A2000055.LST_Day_1km.tif")
# x <- readGDAL("C:/Users/Boris/Documents/MODIS_Windfarm_Temperature/MODIS/PROCESSED/MOD11A1.006_20190312131640/MOD11A1.A2000056.LST_Day_1km.tif")

# bbox_matrix = cbind(c(303090.1, 302182.6, 325722.6, 326522.0, 303090.1),c(6196304, 6239533, 6239995, 6196769, 6196304))
# p = Polygon(bbox_matrix)
# ps = Polygons(list(p),1)
# sps = SpatialPolygons(list(ps))
# proj4string(sps) = CRS("+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# temp = data.frame(f=99.9)
# spdf = SpatialPolygonsDataFrame(sps,temp)
# 
# x = raster(x)
# result = crop(x, spdf)