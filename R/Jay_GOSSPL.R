




# R script file (“GetNearestOceanSoundSpeedProfileForLocation&Month.r”) which provides an example of
# reading from the Rdata sound speed structure and extracting a profile of sounds speed vs. depth for the
# nearest grid point to a given location.

# Load world sound speed data from saved data files and plot & save sound speed profile
# for specific month and nearest locations gridded at 0.25-degree scale
# Sound speed depths are saved for the actual depths in Oceanographic data files and as
# smoothed values at 1-m depth intervals
# Written by Jay Barlow, last modified 8/21/2018
library("gam") #needed to produce smoothed sound speed profile at 1-m scale

# Example location and month
Lat= 33.25; Long= -118.5; iMonth= 7; aName="_ExampleLocation"
setwd("C:/Users/annes/Documents/SoundSpeedWorldwideData")

aMonth= month.abb[iMonth] #abbreviated month names for naming output files
load(paste("WorldSoundSpeedData_",aMonth,".RData",sep=""))
# plot sound speed profiles (south and west are negative lat/long)
iLat= which.min(abs(LatValues-Lat))
iLong= which.min(abs(LongValues-Long))
SoundSpeedMatrix[iLong,iLat,1:78]
imaxDepth= which(is.na(SoundSpeedMatrix[iLong,iLat,1:78]))
imaxDepth= imaxDepth[1]-1
SoundSpeed= SoundSpeedMatrix[iLong,iLat,1:imaxDepth]
Depths= DepthValues[1:imaxDepth]
maxDepth= max(Depths)
plot(SoundSpeed,-Depths,xlab="Sound Speed (m/s)",ylab="Depth",main=(paste("Sound Speed Profile
for",aMonth,"Lat=",LatValues[iLat],"Long=",LongValues[iLong],sep=" ")))
# smooth and plot the sound speed profile for given location
gamout= gam(SoundSpeed~s(Depths,30))
newdata= data.frame(Depths=0:maxDepth)
predicted= as.numeric(predict(gamout,newdata=newdata))
par(cex=1.3)
plot(SoundSpeed,-Depths,pch=19,xlab="Sound Speed (m/s)",ylab="Depth",main=(paste("Sound Speed Profile
for",aMonth,"Lat=",LatValues[iLat],"Long=",LongValues[iLong],sep=" ")))
lines(predicted,-(0:maxDepth), col="red", lwd= 2)
setwd("S:\\SoundSpeedWorldwideData\\OutputProfileFiles")
SoundSpeedByDepth= data.frame(Depth=Depths,SoundSpeed=SoundSpeed)
write.csv(SoundSpeedByDepth,file=paste("SoundSpeedByDepth",aName,".csv",sep=""),row.names=FALSE)
SmoothedSoundSpeedByDepth= data.frame(Depth=0:maxDepth,SoundSpeed=predicted)
write.csv(SmoothedSoundSpeedByDepth,file=paste("SmoothedSoundSpeedByDepth",aName,".csv",sep=""),row.n
          ames=FALSE)
cat("Nearest Gridded Location=",LatValues[iLat],LongValues[iLong])
    