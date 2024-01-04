# Katie's WBM
# remotes::install_github("mikejohnson51/AOI") # suggested!
# remotes::install_github("mikejohnson51/climateR")
library(tidyverse)
library(climateR)
library(AOI)
library(here)
library(plotrix)
library(zoo)
library(ggplot2)
library(grid)
library(cowplot)
library(reshape2)
library(raster)
library(ncdf4)
library(reshape2)
library(WriteXLS)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(forcats)
library(openxlsx); 
library(sf)
library(raster)
library(rgdal)
library(R.utils)
library(tmap); 
library(tmaptools)
library(ggmap)
library(ggspatial)
library(gridExtra)
library(SPEI)
library(sp)
library(skimr)
library(cft)
library(stringr)
library(ggpubr)
library(lemon)
library(ggfortify)
library(extRemes)
library(data.table)
library(sp)
library(terra)

# KW:
# BRCA CURRENT EAST CREEK WELLS SUPPLY
# -112.2124, 37.63281
# -112.2115, 37.63198
# BRCA POTENTIAL FUTURE SUPPLY
# -112.2188, 37.67128


SiteID = "CurrentSupply1"
Park = "BRCA"  
SiteName = SiteID #optional - set at SiteID if don't want to use
Lat = 37.63281
Lon = -112.2124
set.seed(50)

rss_dir <- "RCF_Testing"

local_rss_dir <- paste(rss_dir, Park, sep = "/")

if(dir.exists(local_rss_dir) == FALSE){
  dir.create(local_rss_dir)
}

OutDir <- paste(rss_dir, Park, SiteName, sep = "/")

if(dir.exists(OutDir) == FALSE){
  dir.create(OutDir)
}

DataDir = paste(OutDir, 'input-data', sep = "/") # for .csv's
if(dir.exists(DataDir) == FALSE){
  dir.create(DataDir)
}


### Pulling local parameters from ...?

proj4 <-"+init=epsg:5070" #  North America Albers Equal Area Conic
epsg <- 5070 # North American Albers Equal Area Conic
source('scripts/getParkBoundary.R')

# park <- st_read('data/general/spatial-data/nps_boundary_centroids/nps_boundary_centroids.shp') %>%
#   st_transform(., st_crs(epsg)) %>%
#   filter(UNIT_CODE == Park)

# if 2 units and one has "preserve" in name, do not use
# park <- if(nrow(park) > 1) {
#   park[!grepl("Preserve", park$UNIT_TYPE),]
# } else{park}

# select park
Koppen_parks <- read.csv("data/general/Climate_Zones_by_Unit.csv")

# KW: select state that park falls within... but further down, `s` gets swapped
# KW: to a data fram so removing from the workflow
# s <- st_read('data/general/spatial-data/State_Shapefile/Contig_US_Albers.shp') %>%
#   st_transform(., st_crs(epsg)) %>% filter(STATE_ABBR == Park$STATE)

# TWO DIFFERENT OPTIONS FOR CENTROID - use 1st option if running a general RSS and using park centroid. Second option if using specific lat long.
# state = s$STATE_NAME
# Lat = Koppen_park$Lat
# Lon = Koppen_park$Lon

maca <- raster('data/general/spatial-data/Climate_grid/tdn_90d.nc') %>% 
  projectRaster(., crs = proj4)
dem <- raster('data/general/spatial-data/elevation_cropped.tif') %>%
  projectRaster(., crs = proj4)
soil <- raster('data/general/spatial-data/water_storage.tif') %>%
  projectRaster(., crs = proj4)

nps_boundary <- st_read('data/general/spatial-data/nps_boundary/nps_boundary.shp') %>%
  st_transform(., st_crs(maca))
nps_centroids <- st_read('data/general/spatial-data/nps_boundary_centroids/nps_boundary_centroids.shp') %>%
  st_transform(., st_crs(maca))
US_Counties <- st_read('data/general/spatial-data//US_Counties/tl_2016_us_county.shp') %>%
  st_transform(., st_crs(maca))
US_States <- st_read('data/general/spatial-data/State_Shapefile/Contig_US_Albers.shp') %>%
  st_transform(., st_crs(maca))

# KW: changed place to our SiteID
s <- data.frame(place = SiteID, lat = Lat, lon = Lon)
# s <- data.frame(place = "Water-Supply-1", lat = Lat, lon = Lon)
Site <- st_as_sf(s, coords = c(3,2)) %>% st_set_crs(., 4326)

park <- filter(nps_boundary, UNIT_CODE == Park)
park <- if(nrow(park)>1) {
  park[!grepl("Preserve", park$UNIT_TYPE),]
} else{park}
# in order to use auto zoom feature, must be in lat/long:
state <- US_States %>% filter(STATE_ABBR == park$STATE)

park <- st_transform(park, 4326) 
State <- filter(US_States, STATE_NAME == state$STATE_NAME)
region = park$REGION

state_and_park <- tm_shape(State) +
  tm_borders() + 
  tm_fill(col = "lightgrey") +
  tm_shape(park) +
  tm_borders() + 
  tm_fill(col = "green")

park_and_centroid <- tm_shape(park) + 
  tm_borders() +
  tm_fill(col = "lightgreen") + 
  tm_shape(Site) + 
  tm_dots(size = 1, shape = 19,col="black")

tmap_arrange(state_and_park, park_and_centroid)

# Obtain MACA grid outline (not information within)
Site <- st_transform(Site, st_crs(maca))
park <- st_transform(park, st_crs(maca))
# objects must be Spatial (sp) to work with raster package (cannot be sf):
Site <- as_Spatial(Site)
# find grid cell the location lies within
cell <- cellFromXY(maca, Site) 
maca_cell <- rasterFromCells(maca, cell) # create stand-alone raster for single MACA cell
maca.poly <- rasterToPolygons(maca_cell) # Create MACA polygon - original file in lat/long (note: datum differs from park shapefiles)
maca.poly2<- crop(maca.poly,as_Spatial(park))

tm_shape(park) + 
  tm_borders() + 
  tm_shape(Site) + 
  tm_dots() + 
  tm_shape(maca.poly) + 
  tm_borders()

#####   SLOPE, ASPECT AND RANDOM POINTS   #####

# Create slope and aspect rasters:

slope <- terrain(dem, opt = "slope", unit = "degrees", neighbors = 4) # 4 is better for "smooth" surfaces; 8 is better for rough. See https://www.rdocumentation.org/packages/raster/versions/3.1-5/topics/terrain
aspect <- terrain(dem, opt = "aspect", unit = "degrees")

# get 10 random points from soil raster and create SpatialPoints object
park2 <- maca.poly
# KW: Water supply doesn't always exist in park boundary, so, 
# KW: I'm not cropping cell to park boundary. Original code:
# park.maca.intersect <- crop(park2[1],maca.poly2)
park.maca.intersect <- maca.poly
plot(park.maca.intersect)
points <- spsample(park.maca.intersect, n = 10, type = "random")

# plot to check points appear within borders of MACA cell. 
tm_shape(park) + 
  tm_borders() + 
  tm_shape(maca.poly) + 
  tm_borders() + 
  tm_shape(points) + 
  tm_dots()

####    EXTRACT DATA FROM POINTS  ######################################################################################################

# reproject points to lat/long so can eventually add to .csv

latlong <- st_as_sf(points) %>%
  st_transform(crs = 4326) # project to lat/long

# begin new dataframe for wb_sites:
wb_sites <- as.data.frame(st_coordinates(latlong)) 
wb_sites[,3] <- raster::extract(dem, points)
wb_sites[,4] <- raster::extract(aspect, points)
wb_sites[,5] <- raster::extract(slope, points)
wb_sites[,6] <- raster::extract(soil, points)
wb_sites[,7] <- seq.int(nrow(wb_sites))
wb_sites[,8] <- 5 # default value for wind # KW UPDATE
wb_sites[,9] <- 0 # default value for snowpack # KW UPDATE
wb_sites[,10] <- 0 # default value for Soil.Init # KW UPDATE
wb_sites[,11] <- 1 # default value for shade coefficient # KW UPDATE

wb_sites <- dplyr::select(wb_sites, 7, 2, 1, 3:6, 8:11) # reorder columns
colnames(wb_sites) <- c("WB_site", "Lat", "Lon", "Elev", "Aspect",
                        "Slope", "SWC.Max", "Wind", "Snowpack", "Soil.Init", 
                        "Shade.Coeff")

wb_sites$SWC.Max = wb_sites$SWC.Max*10 # convert units for Soil Water-holding capacity
# KW:
wb_sites$SiteID = SiteID
wb_sites # check to be sure values are populated correctly. There should not be NA values. 

write.csv(wb_sites, file = paste0(OutDir,"/WB_site_parameters.csv"), row.names = FALSE)

# rm(proj4, epsg, dem, soil,US_States,US_Counties,state_and_park,State,slope,points,park_and_centroid,
#    nps_centroids,nps_boundary,latlong,Site,aspect)

CFs_all <- c("Warm Wet", "Hot Wet", "Central", "Warm Dry", "Hot Dry")

##Color schemes

#Colors for CF values plotted side by side (match order of CFs vector)
colors5 <-  c("#6EB2D4", "#05689F", "#F6B294", "#CA0020","grey")
colors5.2 <- c("#6EB2D4", "#05689F", "grey", "#F6B294", "#CA0020")

# Switch for using Tercek csvs or downloading own data:
centroids_csv <- "Y" 
# Switch for method Indiv_method = c("corner", "pca"):
Indiv_method <- "pca" 
# Percentage of models to drop from ranking:
Percent_skill_cutoff = .1 
# Indicates whether Q/I present at bottom of plot to ID CF method used:
MethodCaption = "Y" 

# --------- Information for Climate Futures Analyses (MACA) ------------------ #

# Specify parameters:

Yr = 2050 # Central year
Range = 30  # Number of years to summarize (should be at least 30)
BasePeriod = "1979-2012"
rollLen = 10 # rolling mean length for time series plots

# Threshold percentages for defining Climate futures. Default low/high:  0.25, 0.75
CFLow = 0.25     
CFHigh = 0.75

# Quantiles for temperature threshold calculations
QuantileLow = 0.05
QuantileHigh = 0.95

HotTemp = 95
ColdTemp = 32
PrecipThreshold = 0.05

# ------------ Information for WB analyses ------------------------ #

# Method for PET calculation 
# KW: PET can be massive contributor to projection uncertainty so important for us to explore
# Hamon is default method for daily PRISM and MACA data (containing only Tmax, Tmin, and Date).  KW what about latitude?
# KW: Oudin also only uses air temperature, extraterrestrial radiation, and day length. 

PET_Method <- "Hamon"  

# ------------ Information for drought analyses ------------------------ #

## SPEI variables 
# KW: SPEI = Standardized Precipitation ET Index = drought index

# The value for the period SPEI is aggregated. 6-months is standard but could change if wanted:
SPEI_per <- 6 
# SPEI value under which drought events are considered:
truncation <- -.5 
SPEI_start <- 1980
SPEI_end <- 2012

# Needs aoi to run -- dig back through original cft code at how to create aoi from lat/lon
# AOI <- geocode(location = c("Fort Collins"), pt = TRUE) 
options(timeout = 300)
# AOI <- aoi_get(list(Lat, Lon, .01, .01)) 
# KW: made changes to getting the AOI
AOI <- Site %>%
  st_as_sf() %>%
  st_buffer(., dist = 0.01)

vars = c("tasmax", "tasmin", "pr", "rhsmax", "rhsmin")
# KW a bunch of different climate models:
# This script does not include CCSM or NorESM because they do not
# have Rh and thus break the script. They are both typically central
# tendency models.
gcms = c("bcc-csm1-1",
         "bcc-csm1-1-m",
         "BNU-ESM",
         "CanESM2",
         "CNRM-CM5",
         "CSIRO-Mk3-6-0",
         "GFDL-ESM2G",
         "GFDL-ESM2M",
         "HadGEM2-CC365",
         "HadGEM2-ES365",
         "inmcm4",
         "IPSL-CM5A-LR",
         "IPSL-CM5A-MR",
         "IPSL-CM5B-LR",
         "MIROC5",
         "MIROC-ESM",
         "MIROC-ESM-CHEM",  
         "MRI-CGCM3")

TFtoC <- function(T){(T-32)/1.8}

# VP from FAO -  https://www.fao.org/3/x0490e/x0490e07.htm
# could also use Buck 1981 for 'improved':
# Buck: VPDsat (mb) = (1.0007 + (3.46 * 10^-6 * P)) * 6.1121 * exp((17.50 * T)/(T+240.87))
# where T is deg C and P is atm pressure mb. Above for P > 800 (correction is minimal)
# Zackman re: Ragwala uses: Es = 611.6441 * 10^[(7.591386*T)/(240.7263+T)] where Tavg.
# Shelley sent Vaisala- to use Tavg for 611.6441 parameter - for 020 to +50 C.
# VPsatT = saturation VP @ T deg C [kPa]
# VPD [kPa]
VPsatT <- function(T){0.6108 * exp((17.27 * T)/(T + 237.3))}   

VPD <- function(TminF, TmaxF, RHmin, RHmax){
  Tmin <- TFtoC(TminF); Tmax <- TFtoC(TmaxF)
  es <- (VPsatT(Tmin)+VPsatT(Tmax))/2
  ea <- (VPsatT(Tmin)*RHmax*.01 + VPsatT(Tmax)*RHmin*.01)/2
  es - ea   }  # end VPD  
# scens = c("rcp45", "rcp85")

start.time <-Sys.time()

future_all <- data.frame()

# KW: This took 2 hours for BRCA Well 1 AOI at daily.
# KW: Takes 45 minutes at monthly.
for (i in 1:length(vars)){ 
  for (j in 1:length(gcms)){ 
    print(paste("downloading", gcms[j], vars[5], sep = " "))
    # if(i<3 & j %in% c(5,20)) next
    # cat(i)
    future1 = getMACA(AOI, 
                      #KW changed from daily to monthly per Amber's advice:
                      timeRes = "month",
                      model = gcms[j],
                      varname = vars[i], scenario  = "rcp45", 
                      startDate = "2023-01-01", endDate = "2099-12-31")
    future2 = getMACA(AOI, 
                      #KW changed from daily to monthly per AR's advice:
                      timeRes = "month",
                      model = gcms[j], 
                      varname = vars[i], scenario  = "rcp85",
                      startDate = "2023-01-01", endDate = "2099-12-31")
    if(j == 1) {future <- dplyr::left_join(future1, future2, by = "date")} else{
      try(future <- dplyr::left_join(future, future1, by='date') %>%
        dplyr::left_join(., future2, by='date'))
    }
  }
  
  future_long <- future %>%  
    tidyr::pivot_longer(-date)
  FL <- future_long %>% 
    dplyr::mutate(GCM = gsub("^[^_]*_([^_]+)_.*$", "\\1", future_long$name),
           RCP = sub('.*_', '', future_long$name)) |> 
    dplyr::rename(!!vars[i]:=value) %>% dplyr::select(-c(name))
  if(i==1) { future_all = FL } else {
    future_all = dplyr::left_join(future_all, FL, by=c("date", "GCM", "RCP"))
    # rm(future_long, FL,future1, future2, future)
  }
}

end.time <- Sys.time()

end.time-start.time
# write.csv(future_all,"future_climate_0912.csv",row.names = F)

# future_all <- read.csv("future_ClimateR.csv",header=T)
future_all$date <- as.POSIXct(future_all$date,format="%Y-%m-%d")

Future_all <- future_all %>%
  mutate(Date=date,
         GCM=paste(GCM,RCP,sep="."),
         PrcpIn = pr/25.4,
         TmaxF = ((tasmax-273.15)*9/5) + 32,
         TminF = ((tasmin-273.15)*9/5) + 32,
         RHmaxPct = rhsmax,
         RHminPct = rhsmin,
         TavgF = (TmaxF+TminF)/2,
         Year = format(date,"%Y")) |> 
  dplyr::select(c("Date","GCM","PrcpIn","TmaxF","TminF","RHmaxPct","RHminPct","TavgF","Year", "RCP")) |> 
  mutate(VPD = VPD(TminF, TmaxF, RHminPct, RHmaxPct),
         DOY = yday(Date))

write.csv(Future_all, paste0(OutDir, "/", SiteName,"_future.csv", sep = ''), row.names = F)
Future_all <- read.csv(paste0(OutDir, "/", SiteName,"_future.csv", sep = ''))

####################
#### GRIDMET data
vars = c("tmmx", "tmmn", "pr", "rmax", "rmin")

historical_all <- data.frame()
for (i in 1:length(vars)){
  hist = getGridMET(AOI, varname = vars[i], startDate = "1969-01-01", endDate = "2022-12-31")
  HL <- hist
  if(i==1) { historical_all = HL } else {
    historical_all = left_join(historical_all, HL, by = c("date"))
    rm(HL,hist)
  }
}

Gridmet <- historical_all %>%
  mutate(Date = date,
         GCM = "gridmet.historical",
         PrcpIn = pr/25.4,
         TmaxF = ((tmmx-273.15)*9/5) + 32,
         TminF = ((tmmn-273.15)*9/5) + 32,
         RHmaxPct = rmax,
         RHminPct = rmin,
         TavgF = (TmaxF+TminF)/2,
         Year = format(date,"%Y")) %>% 
  dplyr::select(c("Date","GCM","PrcpIn","TmaxF","TminF","RHmaxPct","RHminPct","TavgF","Year")) |> 
  mutate(VPD = VPD(TminF, TmaxF, RHminPct, RHmaxPct),
         DOY = yday(Date))

write_csv(Gridmet, paste0(OutDir, "/", SiteName,"_historical.csv"))
Gridmet <- read_csv(paste0(OutDir, "/", SiteName,"_historical.csv"))

rm(historical_all, future_all, future_long, FL,future1, future2, future)

# KW: Pull in the plotting functions:
source("scripts/Plotting_Functions.R", local = knitr::knit_global())
# KW: this renaming is weird and I want to change:
SiteID=Park
# KW: develops a lot of summary climate stats for hist and future
# KW: this one doesn't seem to work the correct way...
source("scripts/Plot_Table_Creation.R", local = knitr::knit_global()) # Requires "PARK_init_parsed.RData"; Output = "PARK_lat_long_Final_Environment.RData". 
# KW: this renaming is weird and I want to change:
SiteID=SiteName
# KW: plots of predicted changes in temp and precip across climate futurues:
source("scripts/Scatter_and_diagnostic.R", local = knitr::knit_global()) # Creates scatter and diagnostic plots. Outputs to Figs MACA folder. 

# Specify Climate Futures

## Run this section if already downloaded data objects
# OutDir <- "C:/Users/arunyon/3D Objects/Local-files/Git-repos/CCRP_automated_climate_futures/PARK/NIOB_21/"
DataDir = paste0(OutDir,'/input-data/') # for .csv's
#####################################################
DataFile <- list.files(path = DataDir, pattern = 'Final_Environment.RData', full.names = TRUE) # Environment needs to be added if not parsing MACA data
load(DataFile)

SiteID <- SiteName
FutureSubset <- CFs_all[c(1,5)]; CFs = FutureSubset  # Pick pair of climate futures.
CF_abbreviation <- "WW-HD"
WB_GCMs <- subset(WB_GCMs, CF %in% CFs)

colors2<- colors5[c(1,4)] # Select pair of climate futures - WarmWet/HotDry
#colors2<- c("#F3D3CB","#12045C")  # Select pair of climate futures - HotWet/WarmDry

colors3<-c("white",colors2)
col<- c("darkgray",colors2)  # WarmWet/HotDry
#col<- c("darkgray","#F3D3CB","#12045C")  # HotWet/WarmDry

CFDir = paste0(OutDir,"WarmWet_HotDry/") # for .csv's
if(dir.exists(CFDir) == FALSE){
  dir.create(CFDir)
}

TableDir = paste0(CFDir,"/tables/") # for .csv's
if(dir.exists(TableDir) == FALSE){
  dir.create(TableDir)
}

FigDir = paste0(CFDir,"/figures/") # for .csv's
if(dir.exists(FigDir) == FALSE){
  dir.create(FigDir)
}

if (exists("FutureSubset") == FALSE) stop("Please specify Climate Futures") # At the top of the script, please enter a combination of futures, e.g. Warm Wet/Hot Dry

source(here::here("scripts", "Threshold_Bar_Charts.R")) # Requires "PARK_lat_long_Final_Environment.RData". Outputs plots and Excel Workbook 

# KW: Using climate data to predict additional water-related characteristics (with NPS WBM)
# source(here::here("scripts","Daily_WB_Model.R"))
# 
# source(here::here("scripts", "Drought_characteristics.R"))

# source(here::here("scripts", "WBgrid_plotting.R")) # Run the Water Balance Model

#source(here::here("scripts", "Return_Events.R"))

# source(here::here("scripts", "Summary plots.R"))

# source(here::here("scripts", "Report_plots.R"))
