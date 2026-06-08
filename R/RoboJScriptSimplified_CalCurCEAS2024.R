# Full simple workflow
#### PAMpal processing ####
library(PAMpal)
library(here)

db <- 'F:/Beaked_whale_data/Databases'
bin <- 'F:/Beaked_whale_data/Binaries'
pps <- PAMpalSettings(db=db, binaries=bin, sr_hz='auto', winLen_sec=.0025,
                      filterfrom_khz=10, filterto_khz=NULL)

# data <- processPgDetections(pps, mode='db', id='RoboJStudy')
# data <- setSpecies(data, method='pamguard')
# saveRDS(data,file='F:/Beaked_whale_data/AcousticStudy/AcousticStudy_01_EventsOnly.rds')

# # save this "data" object with saveRDS so you can load later and skip this portion
# recFolder <- 'D:/'
# data <- addRecordings(data, folder=recFolder)
# saveRDS(data,file='F:/Beaked_whale_data/AcousticStudy/AcousticStudy_02_EventswRecordings.rds')
# 
# data <- addGps(data) #GPS tracks in CalCurCEAS 2024 Drifter Analysis Github repo
# saveRDS(data,file='F:/Beaked_whale_data/AcousticStudy/AcousticStudy_03_EventswRecordingsGPS.rds')

#Load in saved Acoustic Study with GPS data
data<-readRDS('F:/Beaked_whale_data/AcousticStudy/AcousticStudy_03_EventswRecordingsGPS.rds')

#### parameters for rest of analysis ####
hpDepth <- 150          # hydrophone depth (m)
animalDepth <- 1217     # mean animal depth (m)
animalSd <- 354         # sd animal depth (m)
truncDist <- 4000       # truncation distance (m)
snapshot <- 2           # snapshot length (minutes)
sdAngleErr <- .59       # sd received angle errors (degrees)
maxTime <- 62           # max dive time (minutes)
ddpMean <- 191.4        # dive depth period mean (minutes)
ddpSd <- 28.2           # ddp sd (minutes)
gsizeMean <- 1.9        # group size mean
gsizeSd <- .13          # group size sd
studyArea <- 1057925                    # study area size. Used to convert density to abundance (km^2)
stationPattern <- '.*_0*([1-9][0-9]*)\\.sqlite3' # pattern to turn DB name into station name
# stationPattern <- '^(DS[0-9]{1,2}).*'  
outPath <- 'RoboJayOutputs'             # path for storing some outputs

#### data processing for RoboJ ####
source('R/RoboJayFunctions.R')
eventClicks <- formatClickDf(data, stationPattern = stationPattern)
eventClicks$station <- as.integer(eventClicks$station)
eventClicks <- splitSSPGroups(eventClicks, splitBy='station', time=3600*24*1)

# # This creates a SSP for each row, has to download from HYCOM server so may be slow
# # Will save the downloaded data to "file" so that it can be read back in and used in the future
# sspList <- createSSPList(eventClicks,
#                          file=file.path(outPath, paste0(as.character(Sys.Date()), '_sspList.rds')),
#                          timeout=360)

# Can read this from disk if running again in future on same data (change date for future)
sspList <- readRDS(here('RoboJayOutputs','2026-06-05_sspList.rds'))

# filter to whatever species for this project
eventClicks <- filter(eventClicks, species == 'ZC')

# Uses the SSP to correct angles using raytrace, then applies raytrace angle correction
# to existing angles. May be quite slow since raytracing can take some time
eventClicks <- doAngleCorrection(eventClicks, sspList=sspList,hpDepth=hpDepth, animalDepth=animalDepth)

#Omit corrected angles that are NA values
eventClicks <- eventClicks[!is.na(eventClicks$angleCorrected), ]
table(is.na(eventClicks$angleCorrected))  # confirm clean

# make this for your dataset
# needs these columns, "station" must match labels created with "stationPattern" earlier
# dutyCycle is "minutes on/minutes off" so 2 of 10 minutes is 2/8
recorderInfo <- tribble(
  ~station, ~dutyCycle, ~recorder, ~deployTime, ~retrTime,
  1, '6/6', 'ST640', '2024-08-17 18:58', '2024-08-25 13:30',
  2, '6/6', 'ST640', '2024-08-20 13:32', '2024-08-30 1:50',
  3, '6/6', 'ST640', '2024-08-20 20:35', '2024-08-30 5:30',
  4, '6/6', 'ST640', '2024-08-21 3:56', '2024-08-28 3:15',
  6, '6/6', 'ST640', '2024-08-23 0:46', '2024-08-30 15:39',
  8, '6/6', 'ST640', '2024-09-13 22:03', '2024-09-21 23:40',
  9, '6/6', 'ST640', '2024-09-14 13:54', '2024-09-23 21:33',
  10, '6/6', 'ST640', '2024-09-14 21:18', '2024-09-24 1:30',
  11, '6/6', 'ST640', '2024-09-15 16:55', '2024-09-22 15:14',
  12, '6/6', 'ST640', '2024-10-01 23:00', '2024-10-05 3:17',
  13, '6/6', 'ST640', '2024-10-02 0:53', '2024-10-04 4:15',
  14, '6/6', 'ST640', '2024-10-02 1:46', '2024-10-05 2:10',
  15, '6/6', 'ST640', '2024-10-05 14:24', '2024-10-11 14:48',
  16, '6/6', 'ST640', '2024-10-06 1:58', '2024-10-13 16:09',
  17, '6/6', 'ST640', '2024-10-06 14:38', '2024-10-12 14:15',
  18, '6/6', 'ST640', '2024-10-06 23:46', '2024-10-11 21:51',
  19, '6/6', 'ST640', '2024-10-09 4:53', '2024-10-15 21:55',
  20, '6/6', 'ST640', '2024-10-26 14:49', '2024-11-09 9:36',
  21, '6/6', 'ST640', '2024-10-27 1:03', '2024-11-07 14:55',
  22, '6/6', 'ST640', '2024-10-27 14:08', '2024-11-09 7:38',
  23, '6/6', 'ST640', '2024-10-28 0:00', '2024-11-11 6:24',
  24, '6/6', 'ST640', '2024-11-11 1:55', '2024-11-24 15:19',
  25, '6/6', 'ST640', '2024-11-20 16:41', '2024-12-03 17:17',
  26, '6/6', 'ST640', '2024-11-22 14:34', '2024-12-02 23:09',
  27, '6/6', 'ST640', '2024-11-23 15:59', '2024-12-03 6:25'
)
recorderInfo$deployTime <- as.POSIXct(recorderInfo$deployTime, format="%Y-%m-%d %H:%M", tz='UTC')
recorderInfo$retrTime <- as.POSIXct(recorderInfo$retrTime, format="%Y-%m-%d %H:%M", tz='UTC')

eventSummary <- formatEventSummary(eventClicks, recorderInfo = recorderInfo, snapshot=snapshot)

snapCount <- tallySnapshots(data, stationPattern=stationPattern, recorderInfo=recorderInfo, snapshot=snapshot)
# can test with doJackknife=FALSE and nSim=1e6, but full run should be TRUE and 1e7
detFunction <- newEstDetFunction(eventSummary, subsetCol = 'recorder',
                                 doJackknife = TRUE, jk_nSamples = 20, nSim=1e6, progress=TRUE,
                                 meanDepth = animalDepth,sdDepth = animalSd, truncDist = truncDist,
                                 sdAngleErr = sdAngleErr, hpDepth = hpDepth)
# plots estimated det function to see if it makes sense
plotDetFun(detFunction, truncDist=truncDist)

detHist <- createDetectionHistory(eventSummary, snapshot = snapshot, maxTime=maxTime,meanDepth=animalDepth,
                                  sdDepth = animalSd, hpDepth = hpDepth, truncDist = truncDist, plot=FALSE)

#### Bayesian model steps ####
source('R/RoboJeffFunctions.R')
ehData <- ehDataPrep(detHist, plot=F, plotDir = file.path(outDir, 'EHPlot'),meanDepth = animalDepth,
                     hpDepth = hpDepth, truncDist = truncDist, maxTime = maxTime)
# currently using open BUGS because that's what Jeff used, this could easily
# be done in whatever flavor of bayes you like
library(R2OpenBUGS)

ehSpecs <- list(
    ni = 10000,  # chain length, including burn-in
    nt = 2,      # thinning rate
    nb = 2500,   # burn-in (e.g., 4000 * thin-rate 5 = 20K samples per chain)
    nc = 2       # number of chains
)

##### Set parameters to monitor #####
ehParams <- c("phi1", "p1",
             "mean.availtime.min")

##### Source model file #####
ehModel <- makeEhModel()

##### Run BUGS analysis #####
# With R2OpenBUGS this opens up the openbugs program while running, you
# will need to exit from there when it is done running
ehOut <- bugs(data=ehData, inits=NULL,
              parameters.to.save=ehParams,
              model.file=ehModel, n.chains=ehSpecs$nc, n.iter=ehSpecs$ni,
              n.burnin=ehSpecs$nb, n.thin=ehSpecs$nt,
              saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL, working.directory=getwd(), clearWD=TRUE)

erSpecs <- list(
    ni = 200000,  # chain length, including burn-in
    nt = 5,     # thinning rate
    nb = 50000,   # burn-in
    nc = 2       # number of chains
)

##### source the model file #####
erModel <- makeErModel()

##### Set parameters to monitor #####
erParams <- c("v.mean", "er.mean",                                 # effective search area and radius
              "logDgrp.hyper", "logDgrp.sig",                      # random effect parameters for log-density
              "logDgrp", "logDgrp.err", "Dper1000", "mu",          # parameters for individual DASBRs
              "mean.of.Errs", "SD.of.Errs", "mean.of.logDgrp",     # descriptive summary stats for individuals DASBR estimates (model checks)
              "Dgrp.mean", "D.mean.per1000", "N",                  # mean density and abundance for the study area (use)
              "Dgrp.mean.II","D.mean.per1000.II",  "N.II",         # mean density and abundance for the study area, alternate (do not use)
              "s", "g0")

erData <- erDataPrep(eventSummary, snaps=snapCount, snapshot=snapshot,
                     detFunction = detFunction$All,
                     ehOut=ehOut,
                     ddpMean=ddpMean,
                     ddpSd=ddpSd,
                     gsizeMean=gsizeMean,
                     gsizeSd = gsizeSd,
                     plot=FALSE, plotDir=file.path(outDir, 'ERPlot'),meanDepth = animalDepth,
                     hpDepth = hpDepth, truncDist = truncDist, A=studyArea)

##### Run BUGS analysis #####
erOut <- bugs(data=erData, inits=NULL, parameters.to.save=erParams, model.file=erModel,
              n.chains=erSpecs$nc, n.iter=erSpecs$ni, n.burnin=erSpecs$nb,
              n.thin=erSpecs$nt, saveExec=FALSE, debug=TRUE, OpenBUGS.pgm=NULL,
              working.directory=getwd(), clearWD=TRUE, bugs.seed=2)
