#Analysis of Parturition Indicators (API)

##"Evaluating alternate methods for estimating the timing of parturition in mule deer"
##By Tabitha A. Hughes, Randy T. Larsen, Kent R. Hersey, Madelon van de Kerk, and Brock R. McMillan##
#created 2021-04-21
#Updated 2024-1-11
##If you use this code, please cite our paper##

#########load required packages#########
library(caTools)
library(adehabitatHR)
library(runner)
library(sp)
library(rgdal)
library(dplyr)


#####read in data and prepare to loop through each individual##########
GPS_data = read.csv("MD_ExampleData.csv", header = T)
MDlist = split(GPS_data, GPS_data$uniqueID)
DeerID = as.character(unique(GPS_data$uniqueID))
resultsN = list()

#####Establish parameter thresholds#######
#parameters may be adjusted according to the average velocity, turning angle, and home range size exhibited immediately after parturition by a training population#
k = 36 #k refers to the window size (number of datapoints used in rolling calculations
velocity = 0.05
turning = 1.8
rMCP = 30

########################run loop################################
####################################################################
for(i in DeerID){
  Deer = MDlist[[i]]
  #change lat/long to UTM coordinates:
  LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
  }
  UTM=LongLatToUTM(Deer$longitude, Deer$latitude,12)
  Deer$X=UTM$X
  Deer$Y=UTM$Y
  #format the date into a format r can use:
  #make sure that the date format matches that of your data
  z= as.POSIXct(strptime(Deer$DateTime,format="%m/%d/%Y %H:%M"))
  Deer$z=z
  #sort by date:
  Deer = Deer[order(Deer$uniqueID, z),]
  
  #there could be duplicates, so remove those if present with the code below:
  Deer$unique=paste(Deer$uniqueID, Deer$z, sep=" ")
  Deer=Deer[order(Deer$unique),]
  Deer=Deer[!duplicated(Deer$unique),]
  
  #create trajectories:#
  library(adehabitatLT)
  t <- as.ltraj(xy = Deer[,c("X","Y")], date = Deer$z, id = Deer$uniqueID)
  #unpack into dataframe:
  d=ld(t) 
  ### The folloing code ensures only data within a two hour interval is included
  d$hour <- as.integer(format(d$date, "%H"))
  d$minute <- as.integer(format(d$date, "%M"))
  MDreg <- d %>%
    filter(minute == 0)
  condition <- ifelse(MDreg$dt != 7200 & MDreg$hour %% 2 != 0, TRUE, FALSE)
  # Use the filter function to remove rows based on the condition
  MDreg <- MDreg %>%
    filter(!condition)
  row.names(MDreg)<-NULL
  MDreg <- na.omit(MDreg)
    
    #######calc rolling Velocity#########
    MDreg$dist = MDreg$dist / 1000
    MDreg$dt = MDreg$dt / 3600
    MDreg$V = MDreg$dist/MDreg$dt
    MDreg$RV = runmean(MDreg$V, k=k, align="left")
    MDreg$julian = strftime(MDreg$date, format = "%j")
    
    #########calc rolling  turning angle######
    MDreg$rel.angle = abs(MDreg$rel.angle)
    MDreg$turnangle = runmean(MDreg$rel.angle, k = k, align = "left")
    
    #######calc Minimum convex polygons (MCP)############
    library(sp)
    MCP = runner(
      x = MDreg,
      k = k,
      f <- function(x){
          x[!is.na(x$x) & !is.na(x$y),]
          x[, c("id", "x", "y")]
          coordinates(x) <- c("x", "y")
          proj4string(x) <- CRS( "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")
          mcp.area(x, percent = 95)
        },
      na_pad = TRUE
    ) 
    MDreg$MCP = NA
    MDreg$MCP = matrix(unlist(MCP), nrow = length(MCP), byrow = TRUE)
    ####shift rMCP values up##########
    shift <- function(x, n){
      c(x[-(seq(n))], rep(NA, n))
    }
    MDreg$MCP <- shift(MDreg$MCP, k)
    ####remove buffers#####
    MDreg = tail(MDreg, -k)
    MDreg = head(MDreg, -k)
    library(dplyr)
    a = filter(MDreg, MDreg$RV < velocity)
    a = filter(a, a$turnangle > turning)
    a = filter(a, a$MCP < rMCP)
    b = a$date[1]
    b = format(b, "%j")
    if(is.na(b)){
      resultsN[[i]] = "Not pregnant"; print(paste("Parturition status identified for", unique(Deer$uniqueID)))
    } else {
      resultsN[[i]] = as.numeric(b); print(paste("Parturition date identified for", unique(Deer$uniqueID)))
    }
    print(b) #returns julian day of parturition
}


