# Full simple workflow
#### PAMpal processing ####
library(PAMpal)
db <- 'path/to/dbs'
bin <- 'path/to/binaries'
pps <- PAMpalSettings(db=db, binaries=bin, sr_hz='auto', winLen_sec=.0025,
                      filterfrom_khz=10, filterto_khz=NULL)
data <- processPgDetections(pps, mode='db', id='RoboJStudy')
data <- setSpecies(data, method='pamguard')
data <- addGps(data)
recFolder <- 'path/to/recordings'
data <- addRecordings(data, folder=recFolder)
# save this "data" object with saveRDS so you can load later and skip this portion

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
stationPattern <- '^(DS[0-9]{1,2}).*'   # pattern to turn DB name into station name
outPath <- 'RoboJayOutputs'             # path for storing some outputs

#### data processing for RoboJ ####
source('R/RoboJayFunctions.R')
eventClicks <- formatClickDf(data, stationPattern = stationPattern)
eventClicks <- splitSSPGroups(eventClicks, splitBy='station', time=3600*24*1)

# This creates a SSP for each row, has to download from HYCOM server so may be slow
# Will save the downloaded data to "file" so that it can be read back in and used in the future
sspList <- createSSPList(eventClicks,
                         file=file.path(outPath, paste0(as.character(Sys.Date()), '_sspList.rds')),
                         timeout=240)

# Can read this from disk if running again in future on same data (change date for future)
# sspList <- readRDS(file.path(outPath, '2022-04-20_sspList.rds'))

# Uses the SSP to correct angles using raytrace, then applies raytrace angle correction
# to existing angles. May be quite slow since raytracing can take some time
eventClicks <- doAngleCorrection(eventClicks, sspList=sspList,hpDepth=hpDepth, animalDepth=animalDepth)

# make this for your dataset
# needs these columns, "station" must match labels created with "stationPattern" earlier
# dutyCycle is "minutes on/minutes off" so 2 of 10 minutes is 2/8
recorderInfo <- tribble(
    ~station, ~dutyCycle, ~recorder, ~deployTime, ~retrTime,
    # 1, '2/8', 'ST4300', '2019-08-05 19:00:00', '2019-08-10 19:00:00',
    # 2, '2/8', 'ST4300', '2019-09-05 19:00:00', '2019-09-12 19:00:00'
)
recorderInfo$deployTime <- as.POSIXct(recorderInfo$deployTime, format='%Y-%m-%d %H:%M:%S', tz='UTC')
recorderInfo$retrTime <- as.POSIXct(recorderInfo$retrTime, format='%Y-%m-%d %H:%M:%S', tz='UTC')

eventSummary <- formatEventSummary(eventClicks, recorderInfo = recorderInfo, snapshot=snapshot)
# filter to whatever species for this project
eventSummary <- filter(eventSummary, species == 'ZC')

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
