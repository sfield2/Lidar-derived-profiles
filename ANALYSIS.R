################################################################################
################################################################################
                       #ANALYSIS FOR FIELD 2023
                      #LIDAR-DERIVED-ROAD-PROFILES
                  #https://doi.org/10.1017/aap.2022.31
################################################################################
################################################################################

# Due to sensitive location information first three steps are built as an example.
# To execute steps 1-3, import elevation raster data set and geospatial transect lines
# Publicly available data developed by author "CHACO_ROAD_PROFILE_DATABASE" and
# CHACO_ROAD_COMPOSITE" can be incorporated starting with step 5 and step 6.


################################################################################
################################################################################
                          # SET UP SCRIPT #
################################################################################
################################################################################
# import libraries
library(raster)
library(rgdal)
library(dplyr)
library(rgeos)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(grid)
library(ggrepel)
library(plyr)
library(tidyr)
library(RColorBrewer)
library(ggbeeswarm)

#Set working directory
setwd("C:/Users/...")


################################################################################
################################################################################
                         #STEP #1: IMPORT DATA#
################################################################################
################################################################################

# Import shapefile line of road
road <- readOGR(dsn="./DATA", layer="ANY_SHAPEFILE_LINE_OF_ROAD")

# Import lidar-derived elevation raster dataset
dem <- raster('./DATA/ANY_DEM.tif',na.rm=T)

################################################################################
################################################################################
            #STEP #2: CREATE TRANSECT LINES#
################################################################################
################################################################################
# Step 2.1: Build transects across road by first creating parallel lines at fixed 
# distance from road

t<-as.data.frame(geom(road))
x <- t$x
y <-t$y
d <- -25   # distance away from the road

segment.shift <- function(x, y, d){
  v <- c(x[2] - x[1],y[2] - y[1])
  v <- v/sqrt((v[1]**2 + v[2]**2))
  vnp <- c( -v[2], v[1] )
  return(list(x =  c( x[1] + d*vnp[1], x[2] + d*vnp[1]), 
              y =  c( y[1] + d*vnp[2], y[2] + d*vnp[2])))
}

xn <- numeric( (length(x) - 1) * 2 )
yn <- numeric( (length(y) - 1) * 2 )

for ( i in 1:(length(x) - 1) ) {
  xs <- c(x[i], x[i+1])
  ys <- c(y[i], y[i+1])
  new.s <- segment.shift( xs, ys, d )
  xn[(i-1)*2+1] <- new.s$x[1] ; xn[(i-1)*2+2] <- new.s$x[2]
  yn[(i-1)*2+1] <- new.s$y[1] ; yn[(i-1)*2+2] <- new.s$y[2]
}

coords <-cbind(xn,yn)

offsetlines <- SpatialLines(list(Lines(list(Line(coords)), ID = 1)))
offsetlines1<-SpatialLinesDataFrame(offsetlines, data.frame(id=1:length(offsetlines)))

# run same analysis with opposite distance
d <- 25   
segment.shift <- function(x, y, d){
  v <- c(x[2] - x[1],y[2] - y[1])
  v <- v/sqrt((v[1]**2 + v[2]**2))
  vnp <- c( -v[2], v[1] )
  return(list(x =  c( x[1] + d*vnp[1], x[2] + d*vnp[1]), 
              y =  c( y[1] + d*vnp[2], y[2] + d*vnp[2])))
}

xn <- numeric( (length(x) - 1) * 2 )
yn <- numeric( (length(y) - 1) * 2 )

for ( i in 1:(length(x) - 1) ) {
  xs <- c(x[i], x[i+1])
  ys <- c(y[i], y[i+1])
  new.s <- segment.shift( xs, ys, d )
  xn[(i-1)*2+1] <- new.s$x[1] ; xn[(i-1)*2+2] <- new.s$x[2]
  yn[(i-1)*2+1] <- new.s$y[1] ; yn[(i-1)*2+2] <- new.s$y[2]
}

coords <-cbind(xn,yn)

offsetlines <- SpatialLines(list(Lines(list(Line(coords)), ID = 1)))
offsetlines2<-SpatialLinesDataFrame(offsetlines, data.frame(id=1:length(offsetlines)))

# Step 2.2: Create series of points along each parallel offset line and draw lines
# between points 
points <- gLength(offsetlines1) # if 1-m spacing of transects are needed use code, if not, 
                                # if not, alter this line by dividing or 
                                #  multiplying gLength
p1<-spsample(testlines1, n = points, type = "regular")
p2<-spsample(testlines2, n = points, type = "regular")

coords <- cbind(p1@coords,p2@coords)

linestrings <- st_sfc(lapply(1:nrow(coords),function(i){
  st_linestring(matrix(coords[i,],ncol=2,byrow=T))
}))

linestrings1 <-as.data.frame(linestrings)
transects <- as_Spatial(linestrings)

transects_final<-SpatialLinesDataFrame(transects, data=linestrings1,match.ID = F)

transects_final$geometry<-1:length(transects)

# Export if you want to inspect transect lines over DEM in other GIS 
# writeOGR(transects_final,"./OUTPUT","transects_final",driver="ESRI Shapefile")


################################################################################
################################################################################
               #STEP #3: EXTRACT ELEVATION VALUES AND NORMALIZE#
################################################################################
################################################################################

# re-import shapefile transects and label with id numbers
lines <- transects_final
lines$id<-1:length(lines)

# Step 3.1: Extract elevation value from every raster cell that is intersected by transects.
# This produces a list of all elevation values intersected for each transect in order
df_l <-raster::extract(dem,lines,along=T)
df<-map_dfr(df_l, ~as_tibble(t(.)))

# Label extracted value list with id numbers so that this data can be merged back with shapefile (data imported as "lines")
# THIS MAY NOT BE NEcESSARY
df$id<-1:length(lines)


# Step 3.2 Normalize transect values. 
# The 2nd-nth elevation values listed in each transect are calculated as a positive or negative change 
# from the first elevation value. The first elevation value is then listed as a "zero" value. 
# This step allows all transects to be directly comparable.

# Create subtraction function which allows all columns to be calculated as positive or negative change from value in first column
minus <- function(first_col,other_col){other_col-first_col}

# create data frame for population of normalized dataset
df_op <- as.data.frame(matrix(ncol=ncol(df),nrow=nrow(df)))

# Run the subtraction function for every row in the extracted value dataset (data created as "df")
# which normalizes all transect outputs
for (i in 1:ncol(df)){
  df_op[,c(i)]<-minus(df[,1],df[,i])
}

# remove artifacts of loop that are not applicable (i.e. final column of output)
df_op<-df_op[1:(length(df_op)-1)]
# list same id for merging back with shapefile
df_op$id<-1:length(lines)
# df_op<-df_op[!is.na(df_op$V1),]

# Step 3.3 Merge normalized output of extracted values with shapefile, and export
finaldf <- merge(lines,df_op,by=c("id"))
writeOGR(finaldf,"./OUTPUT","Lines_with_extracts_EXAMPLE",driver="ESRI Shapefile")


################################################################################
################################################################################
                    #STEP #4: FILTER DATA FOR INDIVIDUAL ROADS#
################################################################################
################################################################################
# You can manually filter data in a GIS to select transects that represent most 
# diagnostic road segments and/or segments that have been ground-truthed. 
# If so, re-import shapefile of manually filtered data.

final <- readOGR(dsn="./OUTPUT", layer="Lines_with_extracts_manual_filter_EXAMPLE")
finaldf<-as.data.frame(final) 
finaldf$transect <- sub("^","transect",1:nrow(finaldf))

# Example of automated filter, which removes less diagnostic transects 
# Step 4.1. Filter to a specific road and select portion of transects that 
# you want to filter (i.e.,V1:V70)
df_op <- finaldf %>%
  subset(Road=="ROAD SAMPLE #1")%>%
  select(transect, V1:V70)

# Step 4.2 build summary of transects (mean and std dev at every point of measure)
stats <- c("Mean", "StDev")
num<-1:ncol(df_op)
n<-ncol(df_op)
summarytbl <- as.data.frame(matrix(NA,n,2, 
                                   dimnames=list(num,stats)))%>%
  mutate(Mean = sapply(df_op[,-(n+1)],mean, na.rm=T))%>%
  mutate(StDev = sapply(df_op[, -(n+1)], sd,na.rm=T))

# Step 4.3 retrieve values that are 1 std dev away from mean 
# at 6th point of measure and subset transect database by values 
b1<- summarytbl[6,1] - (summarytbl[6,2])
t1<- summarytbl[6,1] + (summarytbl[6,2])
filter1 <- subset(finaldf,V5< b1|V5 > t1)
filtered_exampleroad <- subset(df_op, !(transect %in% filter1$transect))


# Step 4.4 export automatically filtered data as csv for easier manipulation
write.csv(filtered_exampleroad,"./all_manualandauto_filtered_EXAMPLE.csv", row.names = FALSE)


################################################################################
################################################################################
              #STEP #5: DEVELOP COMPOSITE#
################################################################################
################################################################################
# After manually and systematically filtering transects, develop composite of remaining
# ground-truthed sections

p_summary<-ggplot(data=filtered_exampleroad)+
  geom_smooth(aes(measurement_value_distance,measurement_value))
filter_composite <-ggplot_build(p_summary)$data[[1]]

filter_composite_test <- filter_composite[,13:14]

          ##### IMPORTANT #### 
# You can incorporate composite of Chaco roads here
# filter_composite_test<- read.csv("./DATA/CHACO_ROAD_COMPOSITE.csv",header=T,fileEncoding = 'UTF-8-BOM')

# first calculate slope of composite
filter_composite_slope <- as.data.frame(matrix(nrow=79,ncol=2))
for(i in 1:79){
  a <- filter_composite_test[i+1,2]
  b <- filter_composite_test[i,2]
  c <- filter_composite_test[i+1,1]
  d <- filter_composite_test[i,1]
  filter_composite_slope[i,1] <- ((a-b)/(c-d))
}

filter_composite_slope$V2 <- filter_composite_test[c(2:80),1]



################################################################################
################################################################################
                 #STEP #6: MEASURE SIMILARITY TO COMPOSITE#
################################################################################
################################################################################

##### IMPORTANT #### 
# You can incorporate Chaco Road database here and test against Chaco composite
# finaldf<- read.csv("./DATA/CHACO_ROAD_PROFILE_DATABSE_1.csv",header=T,fileEncoding = 'UTF-8-BOM')


# step 6.1 subset filtered data to two roads being compared
# example pueblo pintado road compared with composite
pp<-subset(finaldf, Road == "Pueblo Pintado Road")


######  compare north road #####
results_n <- as.data.frame(matrix(nrow=nrow(n),ncol=3))
colnames(results_n) <- c("D-Stat","P-Value","transect")


for(i in 1:nrow(n)){
  # sample consecutive series of 10 transects from road
  test<- as.data.frame(n[i,])
  # standardize by length for direct comparison when plotting
  test<-test %>%
    gather(extractionline,measurement_value,V1:V73,factor_key = T)%>%
    select(extractionline,measurement_value,transect,Road,Confidence,No_of_measures)%>%
    mutate(pointofmeasure = rep(1:73,each=nrow(test)))%>%
    mutate(measurement_value_distance = pointofmeasure*(50/No_of_measures))%>%
    drop_na(measurement_value)
  
  # gather summary data of tested transect 
  p.summary<-ggplot(test)+
    geom_smooth(aes(measurement_value_distance,measurement_value),span=0.5)
  p.summary.build<-ggplot_build(p.summary)
  all_summary <- as.data.frame(p.summary.build$data)
  
  #measure slope of tested transects
  test_slope <- as.data.frame(matrix(nrow=79,ncol=2))
  for(j in 1:79){
    a <- all_summary[j+1,2]
    b <- all_summary[j,2]
    c <- all_summary[j+1,1]
    d <- all_summary[j,1]
    test_slope[j,1] <- ((a-b)/(c-d))
  }
  
  test_slope$V2 <- all_summary[2:80,1]
  
  #limit profile to just road surfaces
  test_slope_lim <- test_slope[19:62,]
  composite_slope_lim <- filter_composite_slope[19:62,]
  
  #compare with ks test
  kstest<-ks.test(test_slope_lim$V1,composite_slope_lim$V1)
  results_n[i,1] <-kstest$statistic
  results_n[i,2] <-kstest$p.value
  results_n[i,3] <-test[1,3]
}

################################################################################
################################################################################
                #  STEP 7.1: Determine significant results #
################################################################################
################################################################################
# https://sparky.rice.edu//astr360/kstest.pdf
# c(alpha)sqrt(n1+n2/n1*n2), n = number of samples
# alpha= 0.001 equates to c(alpha) = 1.95

cv <- 1.95*sqrt((44+44)/(44*44))

# determine   amount with D-stats less than critical value and p-value less than 0.05
result_pp_sig <- subset(results_pp,`D-Stat` <= cv & `P-Value` <= 0.05)


################################################################################
################################################################################
               #  STEP 8: merge back with transect shapefiles #
################################################################################
################################################################################
finaldf_withks <- merge(finaldf,results_all_ks,by="transect")
finaldf_withks <- finaldf_withks[,c(1,2,80:81)]
finaldf_with_ks_shp <- merge(final,finaldf_withks,by="id")
writeOGR(finaldf_with_ks_shp,"./OUTPUT","finaldf_with_ks_shp",driver="ESRI Shapefile")


