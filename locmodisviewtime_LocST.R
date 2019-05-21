
#' Calculate local viewtime from MODIS viewtime raster 
#' @description Calculate local view time from MODIS viewtime raster for specific research area, works for MOD_L2 products and MOD11_A1 (as well as MYD versions) 
#' @param fnam MODIS hdf filename, i.e. without filepath
#' @param ra extent of the research area in projection of the research area  (=parameter newproj)
#' @param ltz local time zone, e.g. "Pacific/Auckland" or "EST" to be found in R function OlsonNames()
#' @param viewtime Day_view_time or Night_view_time raster inside .hdf in latlon
#' @param newproj projection used for data in research area 
#' @return data frame with viewtime raster pixel coordinates (lon, lat), pixel value of viewtime raster (lst), 
#' UTC date extracted from the filename (fnamdate), days to be added in case of surpassing 24h (daych), 
#' pixelvalue converted to hour format (lst_24h), the complete date in local solar time (lst_date), the complete 
#' UTC date based on reconversion from local solar time (utc_date) and  the date in the local time zone format converted
#' from the local solar time via UTC (ltz_date) 
#' @author Maite Lezama Valdes
#' @seealso \code{\link{getmodisdate}}
#' @examples
#' ra <- c(390015.8, 392069.4, -1290438, -1288308)
#' data(dvt)
#' viewtime <- dvt
#' fnam <- "MOD11A1.A2018019.h21v16.006.2018020085038.hdf"
#' ltz <- "Pacific/Auckland"
#' newproj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#' locmodisviewtime(fnam, viewtime, ra, ltz, newproj)
#' @export locmodisviewtime
#' @aliases locmodisviewtime






library(raster)


dec_time <- function(dectime){
  h <- sapply(seq(dectime), function(i){
    if(dectime[i] < 1){ # if <1h
      paste0("00:",dectime[i]*60)
    } else if((dectime[i] - floor(dectime[i]))==0){
      paste0(dectime[i], ":00")
    } else {
      paste0(floor(dectime[i]), ":", round((dectime[i]-floor(dectime[i]))*60, digits=0))
    }
  })
  return(h)
}


locmodisviewtime <- function(fnam, viewtime, ra, ltz, newproj){
  
  library(lubridate)
  
  # read latlon viewtime raster 
  vt <- viewtime
  
  # crop viewtime raster to research area 
  e <- extent(ra)
  e <- as(e, "SpatialPolygons")
  proj4string(e) <- newproj
  e.geo <- sp::spTransform(e, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  vt_s <- crop(vt, e.geo)
  vt_s <- vt_s*0.1 # apply scale factor (see 3.2. Scientific Data Sets (SDS) 
  #in https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf)
  
  # get UTC date from fnam
  su <- strsplit(fnam, "A")
  su <- su[[1]][length(su[[1]])]
  org <- paste0(substring(su, 1,4), "-01-01")
  utcday <- as.Date((as.numeric(substring(su, 5,7))-1), origin=org)
  
  if(grepl("v", su)){
    utctime <- paste0(substring(su, 10, 11), ":", substring(su, 13, 14))
  } else {
    utctime <- paste0(substring(su, 9, 10), ":", substring(su, 11, 12))
  }
  
  utcdate <- strptime(paste0(utcday,utctime), format='%Y-%m-%d %H:%M', tz="UTC")
  
  # get location and time values
  xy <- as.data.frame(rasterToPoints(vt_s))
  lst <- xy[,3]
  names(xy) <- c("lon", "lat", "LST")
  
  xy$fnamdate <- utcdate
  
  #+24 if local solar time < 0 or -24 if local solar time >= 24
  
  # my idea about day change: if LST raster contains 20 and 25.5, doesn't it make 
  # sense to assume, that although we substract 24h to get to 1:30, it is the next 
  # day, and not the day before? which makes measurements closer to each other
  
  utc_day <- lapply(seq(lst), function(i){
    if(lst[i] < 0){
      data.frame(daych = -1, lstcorr = lst[i]+24)
    } else if(lst[i] >= 24){
      data.frame(daych = +1, lstcorr = lst[i]-24)
    } else {
      data.frame(daych=0, lstcorr=lst[i])
    }
  })
  
  # daych shows how to change LST day
  utc_d <-do.call(rbind, utc_day)
  xy <- cbind(xy, utc_d)
  
  # make date format for LST
  lst_h <- dec_time(xy$lstcorr)
  xy$lst_h <- lst_h
  
  # convert UTC date from filename to local time zone and take the day from there
  locdate <- with_tz(utcdate, tzone = ltz)
  locday <- format(locdate, "%Y-%m-%d")
  
  # assignment of local timezone is just to evade getting a completely false one in,
  # actually, this is local solar time
  locdates <- strptime(paste0(locday, "_", xy$lst_h), format='%Y-%m-%d_%H:%M', tz=ltz)
  locdates$mday <- locdates$mday+xy$daych
  
  xy$lst_date <- locdates
  
  # Grid is in local solar time, which is: 
  # UTC time plus grid's longitude in degrees / 15 degrees 
  # (in hours, +24 if local solar time < 0 or -24 if local solar time >= 24).
  
  # convert to UTC
  utcconv <- locdates
  # change utcconv to UTC time zone format
  attr(utcconv, "tzone") <- "UTC"
  utcconv$zone <- "UTC"
  
  # find lon / 15 conversion: hours and minutes
  lon_15 <- xy$lon / 15
  lon_15_h <- dec_time(lon_15)
  hrs <- as.numeric(substring(lon_15_h, 1,2))
  mins <- as.numeric(substring(lon_15_h, 4,5))
  
  # convert LST to UTC via latlon/15
  utcconv$hour <- utcconv$hour - hrs
  utcconv$min <- utcconv$min - mins
  xy$utc_from_lst <- utcconv
  
  # change UTC to local time zone
  local_tz_date <- with_tz(xy$utc_from_lst, tzone = ltz)
  xy$ltz_date <- local_tz_date
  
  names(xy) <- c("lon", "lat", "lst", "fnamdate", "daych", "lst_24", "lst_24h", "lst_date", "utc_date", "ltz_date")
  return(xy)
}


# #main <- "E:/MODIS/flight9/MYD11_L2/10_40/"
# main <- "E:/MODIS/flight9/MYD11_L2/"
# 
# list.files(main)
# 
# #fnam <- list.files(main, pattern=".hdf", full.names = F)
# fnam <- "MYD11_L2.A2018019.0905.006.2018020174548.hdf"
# ra <- c(390015.8, 392069.4, -1290438, -1288308)
# viewtime <- raster(list.files(main, pattern="viewtime.tif$", full.names=T))
# #viewtime <- raster(list.files(main, pattern="view_time.tif$", full.names=T))
# ltz <- "Pacific/Auckland"
# newproj <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# 
# 
# locmodisviewtime(fnam, viewtime, ra, ltz, newproj)
# 
