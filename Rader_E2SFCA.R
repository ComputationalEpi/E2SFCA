remove(list = ls())

#load libraries
library(tidyverse)
library(data.table)
library(rgdal)
library(sf)
library(sp)
library(raster)
library(velox)
library(foreach)
library(doParallel)

#set working directory
setwd("C:/Users/benbi/Documents/BU/Boston Childrens/COVID/VF_COVID")

#todays date
date1<-Sys.Date()

#####Read in data####
df3 <- fread("data/va_example_data.csv")

#read in state borders
#source: https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip
statebords<-readOGR("C:/Users/benbi/Documents/BU/Boston Childrens/Vaccine Finder/Vaccine Access/Shapefiles/cb_2018_us_state_500k.shp")

#remove boundaries of areas that we are not studying
statebords2 <- subset(statebords, NAME !="Alaska" & NAME != "Hawaii" & NAME != "Puerto Rico" 
                      & NAME != "Commonwealth of the Northern Mariana Islands"
                      & NAME != "Guam" & NAME != "United States Virgin Islands"
                      & NAME != "American Samoa")

#ensure state data has correct CRS
crs(statebords2) <- "+proj=longlat +datum=NAD83 +ellps=GRS80"

#add convert text lat and long to spatial information for doses
coordinates(df3) <- ~ longitude + latitude

#set common projection
proj4string(df3) <- proj4string(statebords2)

#check
plot(statebords2)
plot(df3, col="red", add = TRUE)

#crop the provider locations by the states we are studying
df4 <- raster::crop(df3, statebords2)

#check to ensure points we exclude are red
plot(df4, col = "green", add = TRUE)

#download the accessibility surface
#source: https://malariaatlas.org/research-project/accessibility-to-cities/
usa <- raster("D:/Friction/2015_friction_surface_v1.geotiff")
usa1 <- crop(usa,statebords2)

#check
plot(usa1)

#read in population count raster
#source https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-count-rev11
worldpop <-raster("C:/Users/benbi/Documents/BU/Boston Childrens/Vaccine Finder/Vaccine Access/grids/gpw_v4_population_count_rev11_2020_30_sec.tif")

#clip to usa
usapop <- crop(worldpop,statebords2)

#mask to not include border regions
usapop2 <- mask(usapop,statebords2)

#check
plot(usapop2)

#convert to velox for quicker count extraction
usapop3<-velox(usapop2)

#set up parallel environment
n <- 2 #number of cores
cl <- makeCluster(n)
registerDoParallel(cl)

####calculations are based on miles and miles per hour. may need to be simplified from km/h
#max travel time you want in minutes
#higher times are more computationally intensive because we caclulate accessability time for a radius of mtd cells (you can move mtd distance in each direction from center)
tt <- 90
#for 90 minutes an invidual can travel x miles, where x is a function of max land speed as allowed by the friction map
maxrs <- 80 #from osm: https://wiki.openstreetmap.org/wiki/OSM_tags_for_routing/Maxspeed

#max travel distance
mtd <- maxrs/60 *tt ##in tt travel minutes, someone can travel at most mtd miles

#convert miles to meters
mrad <- mtd*1609.34

#read in emirical discount function
discount_func <- readRDS("discount_func2.rds")

####Step 1#####
#intialize data and run in parellel
#n is the number of vaccine disribution sites we are estimating for
results1<-foreach(n=1:nrow(df4), .combine=cbind, .packages = c("raster","sf","gdistance","sp","velox")) %dopar%{
point1<-df4[n,] #start with the first site
point2<-buffer(point1, mrad) #draw a radius around it of mrad meters
usa3 <- crop(usa1,point2) #crop the friction map to this radius
T <- gdistance::transition(usa3, function(x) 1/mean(x), 8) #make the transition matrix
T.GC <- gdistance::geoCorrection(T) #geocorrect the transition matrix
access.raster <- gdistance::accCost(T.GC,point1) #run accumulated cost algorithm across the entire friction map to the vaccine site
poly1<-rasterToPolygons(access.raster, fun=function(x){x<=tt}, n=8, na.rm=TRUE, digits=8, dissolve=FALSE) #convert all raster squares under tt minutes to a polygon
poly2<-st_as_sf(poly1) #convert that polygon to a simple feature
poly3<-cbind(poly2,population=usapop3$extract(sp=poly2,fun=mean)) #for each cell in that simple feature, extract the population count
Sum1 <- sum(poly3$population*discount_func(poly3$layer), na.rm=TRUE) #for each cell, apply discount function to layer distance (measured by accumulated cost algorithm). multiply the result by the population count. the total is the denominator of step 1
}

#########
df6<-st_as_sf(df4) #convert to simple feature
df7<-df6[1:length(results1),] #ensure all values were measured (df6 and df7 same length, can use this step to test shorter datasets)
df7$pop <- as.data.frame(t(results1))$V1 #add vector of catchment area population counts to dataframe
df7$rval <- df7$Doses/df7$pop #calculate rval
df7$rval <- ifelse(df7$rval==Inf,0,df7$rval) #if rval is infinite (catchment area population 0), remove influence from dataset

df8 <- df7 %>% dplyr::select(rval) #only look at rval variable
df8 <- st_transform(df8, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") #reproject measure

#2010 census population weighted centroids
tract_weight<-fread("CenPop2010_Mean_TR.txt", 
                    colClasses = c("character","character","character","numeric","numeric","numeric"))

#get full fips from subparts
tract_weight$FIPS <- paste0(tract_weight$STATEFP,tract_weight$COUNTYFP,tract_weight$TRACTCE)

#only look at lower 48 + DC
tract_weight2 <- tract_weight %>% dplyr::filter(!(STATEFP %in% c("02","15","60","66","69","72","78")))

#make spatial
coordinates(tract_weight2) <- ~LONGITUDE+LATITUDE

########## STEP 2
results2<-foreach(n=1:nrow(tract_weight2), .combine=cbind, .packages = c("raster","sf","gdistance","sp","velox")) %dopar%{
  tract1<-tract_weight2[n,] #for each census tract
  tract2<-buffer(tract1, mrad) #draw mrad buffer
  usa3 <- crop(usa1,tract2) #crop the friction map to this radius
  T <- gdistance::transition(usa3, function(x) 1/mean(x), 8) #make the transition matrix
  T.GC <- gdistance::geoCorrection(T)  #geocorrect the transition matrix
  access.raster <- gdistance::accCost(T.GC,tract1)  #run accumulated cost algorithm across the entire friction map to the tract centroid
  poly1<-rasterToPolygons(access.raster, fun=function(x){x<=tt}, n=8, na.rm=TRUE, digits=8, dissolve=FALSE)  #convert all raster squares under tt minutes to a polygon
  poly2<-st_as_sf(poly1) #convert to simple feature
  poly3<-st_join(poly2,df8, join = st_contains, left=TRUE)  #find all vaccine sites within that poly2
  sum2<-sum(poly3$rval*discount_func(poly3$layer), na.rm=TRUE)  #for each cell, apply discount function to layer distance (measured by accumulated cost algorithm). multiply the result by the rval. the total is the accessbiity score for that location
}

#save results!
name2<-paste0("results2_",date1,".csv")
fwrite(results2,name2)