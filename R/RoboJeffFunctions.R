# RoboJeff
library(dplyr)
ehDataPrep <- function(detHist, meanDepth=1217, hpDepth=110, truncDist=4e3, plot=FALSE, plotDir='.', maxTime=61.999) {
    intMatrix <- detHist$dhMat
    intMatrix[intMatrix == -1] <- NA  # if duty-cycle-off was indicated by -1 in the csv, this line changes those to NAs (we want NAs for occasions where detectors were off)
    # Create interval length matrix
    lenMatrix <- matrix(0, nrow=nrow(intMatrix), ncol=ncol(intMatrix))
    for(i in 1:nrow(lenMatrix)) {
        lenMatrix[i, ] <- detHist$dcMap[[detHist$dhDf$station[i]]]
    }
    ##### Create other matrrix objects #####
    # empty matrices of correct dimension
    detMatrix <- effMatrix <- posAngleMatrix <- intMatrix

    # binary detection matrix (changes detection angles to 1s, changes NAs to 0s)
    # if duty-cycle off are NAs in the csv:  matrix that changes angles to 1 (binary indicator of detection) and 0s otherwise
    detMatrix[intMatrix>0] <- 1
    detMatrix[is.na(intMatrix)] <- 0

    # effort matrix (1 if detector on, 0 if off due to duty cycling)
    # (if duty-cycle off are NAs in the csv) matrix of 1s for duty-cyle on, 0s for duty-cycle off
    effMatrix[!is.na(effMatrix)] <- 1
    effMatrix[is.na(effMatrix)] <- 0

    # matrix of positive measure angles (matrix has angles to positive detections and NAs otherwise)
    posAngleMatrix[posAngleMatrix == 0] <- NA

    ##### Derived data and covariate fields #####
    # Estimated distance to initial detection (based on angle to initial encounter), assuming animal is at its mean dive depth
    detDist0_km <- (meanDepth-hpDepth) * tan((180-intMatrix[,1])*pi/180) / 1e3
    # Mean detection distance (average of all detection distances in the encounter)
    detDistMean_km <- (meanDepth-hpDepth) * rowMeans(tan((180-posAngleMatrix)*pi/180),na.rm=T) / 1e3
    # Difference between initial and mean detection
    distDiff <- detDist0_km - detDistMean_km
    # Number of detections per encounter history
    numDet <- rowSums(detMatrix)
    # Number of intervals during which the detector was on (duty cycling), a measure of effort
    numInts <- rowSums(effMatrix)
    # crude measure of detection rate: proportion of intervals detected out of the number of on-intervals
    dpue <- numDet/numInts
    # last interval in which animal was detection (a measure of encounter length)
    lastInt <- apply(detMatrix, 1, function(x) {
        isDet <- x > 0
        if(!any(isDet)) {
            return(0)
        }
        max(which(isDet))
    })
    # alternate crude measure of detection rate: proportion of intervals detected within the period that it was being heard
    dpue2 <- numDet/lastInt

    ##### Create BUGS data variables and object #####
    # only retain records for which initial detection distance is within the truncation distance
    # but note that covariate for analysis is the mean detection distance (for those records where initial distance was < 4km)
    # sample size (number of encounter histories)
    distOK <- detDist0_km < (truncDist/1e3)
    ndives <- nrow(detMatrix[distOK, ])
    # number of encounter history occasions
    nocc <- ncol(intMatrix[distOK, ])
    # index whether DASRR is on 100% of time (dutyidx=1) or duty-cycled (dutyidx=0)
    dutyidx <- rep(NA, ndives)
    for(i in 1:ndives){
        dutyidx[i] = 1 - (sum(effMatrix[distOK,][i, ]==0)>0)
    }
    # BUGS data object
    bugsDataEH <- list(nOcc=nocc,
                       nInd=ndives,
                       eh=detMatrix[distOK, ],
                       dCyc=effMatrix[distOK, ],
                       intLen = lenMatrix[distOK, ],
                       nMin = ceiling(maxTime))
    if(plot) {
        if(!dir.exists(plotDir)) {
            dir.create(plotDir)
        }
        # Plot initial vs. mean detection distances, all data
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km, breaks=50, main=paste("n = ", length(detDist0_km),sep=""))
        hist(detDistMean_km,breaks=50, main=paste("n = ", length(detDist0_km),sep=""))
        plot(detDist0_km,detDistMean_km); abline(a=0,b=1); hist(distDiff, breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to MEAN detection distance < truncDist distance
        meanTrunc <- detDistMean_km < (truncDist/1e3)
        initTrunc <- detDist0_km < (truncDist/1e3)
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[meanTrunc], breaks=25, main=paste("n = ", length(detDist0_km[meanTrunc]),sep=""))
        hist(detDistMean_km[meanTrunc],breaks=25, main=paste("n = ", length(detDistMean_km[meanTrunc]),sep=""))
        plot(detDist0_km[meanTrunc],detDistMean_km[meanTrunc], main="")
        abline(a=0,b=1)
        hist(distDiff[meanTrunc], breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc meanDist 4km.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[initTrunc], breaks=25, main=paste("n = ", length(detDist0_km[initTrunc]),sep=""))
        hist(detDistMean_km[initTrunc],breaks=25, main=paste("n = ", length(detDistMean_km[initTrunc]),sep=""))
        plot(detDist0_km[initTrunc],detDistMean_km[initTrunc], main="")
        abline(a=0,b=1)
        hist(distDiff[initTrunc], breaks=50)
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc initialDist 4km.tiff"), type="tiff")
        dev.off()

        ### Repeat some of the above plots, just for detectors that were never duty cycled

        # Plot initial vs. mean detection distances, data truncated to MEAN detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[meanTrunc & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)",
             main=paste("n = ", length(detDist0_km[meanTrunc & numInts==max(numInts)]),sep=""))
        hist(detDistMean_km[meanTrunc & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)",
             main=paste("n = ", length(detDistMean_km[meanTrunc & numInts==max(numInts)]),sep=""))
        plot(detDist0_km[meanTrunc & numInts==max(numInts)],detDistMean_km[meanTrunc & numInts==max(numInts)],
             xlab="Initial detection km", ylab="Mean detection km", main="")
        abline(a=0,b=1)
        hist(distDiff[meanTrunc & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc meanDist 4km - noDutyCycle.tiff"), type="tiff")
        dev.off()

        # Plot initial vs. mean detection distances, data truncated to INITIAL detection distance < truncDist distance
        x11()
        par(mfrow=c(2,2))
        hist(detDist0_km[initTrunc & numInts==max(numInts) ], breaks=25, xlab="Initial Distance (km)",
             main=paste("n = ", length(detDist0_km[initTrunc & numInts==max(numInts)]),sep=""))
        hist(detDistMean_km[initTrunc & numInts==max(numInts)],breaks=25, xlab="Mean Distance (km)",
             main=paste("n = ", length(detDistMean_km[initTrunc & numInts==max(numInts)]),sep=""))
        plot(detDist0_km[initTrunc & numInts==max(numInts)],detDistMean_km[initTrunc & numInts==max(numInts)],
             xlab="Initial detection km", ylab="Mean detection km", main="")
        abline(a=0,b=1)
        hist(distDiff[initTrunc & numInts==max(numInts)], breaks=50, xlab="Initial distance - Mean distance", main="")
        savePlot(filename = file.path(plotDir, "Detection Distance Summaries - trunc initDist 4km - noDutyCycle.tiff"), type="tiff")
        dev.off()

        ### Plot mean detection distance vs. some encounter rate metrics
        axisSize = 1.5
        labSize = 1.75
        fontType = 2 # bold
        #x11()
        tiff(filename = file.path(plotDir, "Encounter rate metrics vs mean distance.tif"), width=10, height=10,units="in",res=300)
        par(mfrow=c(2,2), mar=c(5,5,4,2))
        # Set 1: all data within truncation range
        plot(detDistMean_km[initTrunc], dpue[initTrunc], main="",
             xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        plot(detDistMean_km[initTrunc], lastInt[initTrunc], main="",
             xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        # Set 2: only for non duty-cycled detectors
        plot(detDistMean_km[initTrunc & numInts==max(numInts)], dpue[initTrunc & numInts==max(numInts)],
             main="", xlab="Mean detection distance (km)", ylab="Detections per on-interval", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        plot(detDistMean_km[initTrunc & numInts==max(numInts)], lastInt[initTrunc & numInts==max(numInts)],
             main="", xlab="Mean detection distance (km)", ylab="Last interval detected", cex.lab=labSize, cex.axis=axisSize, font.lab=fontType)
        mtext(paste("All data within truncation distance (n = ", length(detDistMean_km[initTrunc]), ")", sep=""), side=3, outer=T,cex=2, line= -3, font=2)
        mtext(paste("Only from detectors not duty-cycled (n = ", sum(initTrunc & numInts==max(numInts)), ")", sep=""), side=3, outer=T, line=-33,cex=2,font=2)
        #savePlot(filename = "Encounter rate metrics vs mean distance.tif", type="tiff")
        dev.off()

        # Plot initial detection distance vs. some encounter rate metrics
        x11()
        par(mfrow=c(2,2))
        # Set 1: all data within truncation range
        plot(detDist0_km[initTrunc], dpue[initTrunc], main=paste("n = ", length(detDist0_km[initTrunc]),sep=""),
             xlab="Initial detection distance", ylab="Detections per on-interval")
        plot(detDist0_km[initTrunc], lastInt[initTrunc], main="All data within truncation distance",
             xlab="Initial detection distance", ylab="Last interval detected")
        # Set 2: only for non duty-cycled detectors
        plot(detDist0_km[initTrunc & numInts==max(numInts)], dpue[initTrunc & numInts==max(numInts)],
             main=paste("n = ", sum(initTrunc & numInts==max(numInts)),sep=""), xlab="Initial detection distance", ylab="Detections per on-interval")
        plot(detDist0_km[initTrunc & numInts==max(numInts)], lastInt[initTrunc & numInts==max(numInts)],
             main="Only from detectors not duty-cycled", xlab="Initial detection distance", ylab="Last interval detected")
        savePlot(filename = file.path(plotDir, "Encounter rate metrics vs initial distance.tiff"), type="tiff")
        dev.off()

        # duration length (last interval detected) on non duty-cycled detectors
        mean(lastInt[initTrunc & numInts==max(numInts)])  # mean value for "last Interval", for non-duty cycled detectors, within truncation distance
        # 7.57 for 4k trunction
        length(lastInt[initTrunc & numInts==max(numInts)])  # sample size
        # n=44 for 4k truncation
    }
    bugsDataEH
}

erDataPrep <- function(es, snaps, meanDepth=1217, hpDepth=110, truncDist=4000,
                       esrMean=NULL, esrSd=NULL, detFunction=NULL,
                       A=1057925, ehOut, snapshot=2,
                       ddpMean=191.4, ddpSd=28.2,
                       gsizeMean=1.9, gsizeSd=.13,
                       plot=FALSE, plotDir='ERPlot') {
    truncAngle <- atan(truncDist / (meanDepth-hpDepth)) * 180 / pi
    es <- filter(es, !is.na(meanAngle),
               meanAngle <= truncAngle,
               nClicks >= 3)
    es$station <- as.character(es$station)
    nDets <- group_by(es, station) %>%
        summarise(nDets = n())
    nDets <- left_join(snaps, nDets, by='station')
    nDets$nDets[is.na(nDets$nDets)] <- 0
    nDets <- filter(nDets,
                    !is.na(nSnaps),
                    nSnaps > 0)
    if(plot) {
        if(!dir.exists(plotDir)) {
            dir.create(plotDir)
        }
        # histogram of angles for retained detections
        x11()
        par(mfrow=c(1,2))
        hist(es$meanAngle, breaks=seq(0,90,5), main='Mean Angle')

        # histogram of distances
        hist((meanDepth-hpDepth)*tan(es$meanAngle*pi/180),
             breaks=50, xlab="RANGE (KM)", main="Det Dist if\nAnimal Depth = Mean Depth")
        savePlot(filename = file.path(plotDir, "Detection Angle and Distance Summaries.tiff"), type="tiff")
        dev.off()
        # histogram of detection angles by recorder type
        recs <- unique(es$recorder)
        x11()
        par(mfrow = c(1,length(recs)))
        for(i in seq_along(recs)) {
            hist(es$meanAngle[es$recorder==recs[i]],breaks=seq(0,90,5), xlab="Angle", main=recs[i])
        }
        savePlot(filename = file.path(plotDir, "Detection Angle by Recorder.tiff"), type="tiff")
        dev.off()
    }

    if(is.null(esrMean) &
       !is.null(detFunction)) {
        esrMean <- detFunction$mean_jk_EDR / 1e3
    }
    if(is.null(esrSd) &
       !is.null(detFunction)) {
        esrSd <- esrMean * detFunction$jackKnife_cv
    }
    list(nSites = nrow(nDets),
         n = nDets$nDets,
         A=A,
         k=nDets$nSnaps,
         esrMean=esrMean,
         esrSd=esrSd,
         # See Jay's paper on snapshot length for explanation of adding the duration (Barlow 2021, in press)
         availMean = ehOut$mean$mean.availtime.min + snapshot,
         availSd = ehOut$sd$mean.availtime.min,
         ddpMean=ddpMean,
         ddpSd=ddpSd,
         gsizeMean=gsizeMean,
         gsizeSd=gsizeSd
    )
}

makeEhModel <- function(file=NULL) {
    if(is.null(file)) {
        file <- tempfile()
    }
    cat("

  model{

    for(i in 1:nInd){
      for(t in 2:nOcc){              # n.occ = number of occasions, including initial capture occasion
        phi.it[i,t] <- pow(phi1, (intLen[i,t-1] + intLen[i,t])/2)
        p.it[i,t] <- (1 - pow(1-p1, intLen[i, t])) * dCyc[i,t]     # p = p0 when dcyc[i,t]=1, and equals 0 when dcyc[i,t]=0 (detector is off).  Note, dcyc is a data object.
      }  # t
    }  # i

    ### PRIORS
    phi1 ~ dunif(0,1) # 1 minute
    # phi ~ dunif(0,1) # is estimated per minute in new model
    # p ~ dunif(0,1)
    p1 ~ dunif(0,1) # 1 minute

    ### LIKELIHOOD
    for(i in 1:nInd){
      z[i,1] <- eh[i,1]                           # z indicates whether indiv i is available to sampling; this is 1 in the 1st occasion for all individuals since the encounter history is conditioned on the initial detection
      for(t in 2:nOcc){
        mu.avail[i,t] <- z[i,t-1] * phi.it[i,t]   # mu.avail = prob i survived to t, given it was in sampling population at t-1 (z[i,t-1]=1)
        z[i,t] ~ dbern(mu.avail[i,t])             # index whether animal is still available to detection during capture occasion t
        mu.det[i,t] <- z[i,t] * p.it[i,t]            # mu.det =  prob of being detected (0 if animal no longer avaiable; p otherwise)
        eh[i,t] ~ dbern(mu.det[i,t])              # likelihood for the data (eh[i,t] is the data object)
      } # t
    } # i

    ### DERIVED PARAMETERS

      # Median availability time (in intervals and minutes)
      med.availtime.min <- log(0.5)/log(phi1)     # time in minutes for which 50% of animals are available to detection
      # med.availtime.ints <- med.availtime.min / 2      # number of intervals (inluding 1st) for which 50% of animals are available to detection

      # Mean availablity time (in intervals and minutes)
      # mean.availtime.min <- mean.availtime.ints * 2
      for(i in 1:nMin){
        y[i] <- d[i] * ((i-1) + 0.5)             # y(i) = d(i) * (i + 0.5), where d(i) is the proportion that die between i and i+1
        d[i] <-  pow(phi1, (i-1)) * (1 - phi1)     # d(i) = l(i) * q(i), where l(i) is survivorship to age i; q(i) = 1-phi(i)
      }
      mean.availtime.min <- sum(y[])

  } # end model

", fill=TRUE, file=file)
    file
}

makeErModel <- function(file=NULL) {
    if(is.null(file)) {
        file <- tempfile()
    }
    cat("

    model{

      ##### EFFECTIVE SEARCH RADIUS AND AREA #####
      v.mean <- er.mean * er.mean * 3.14159  # effective search area = pi * r^2
      er.mean  ~ dnorm(esrMean,tau.er)     # r, the effective search radius (3.13 km), from Jay's max sim lik model (for 4km truncation, Jan 2021)
      tau.er <- 1/(sd.er * sd.er)
      sd.er <- esrSd                  # Jay provided CV = 0.10 (so SD = 0.313)

      # This loop accommodates different ESRs for each DASBR, but we ultimately used a single ESR, so er[station] = er.mean
      for(z in 1:nSites){
        er[z] <- er.mean
        v[z] <- 3.14159 * er[z] * er[z]
      } # z


      ##### MODEL FOR DASBR-SPECIFIC DENSITIES #####
      for(z in 1:nSites){                              # loop through DASBR stations
        n[z] ~ dpois(mu[z])                             # liklihood.  Number of *groups* detected per station (n[z] are the data)
        mu[z] <- Dgrp[z] * v[z] * k[z]    # mu = D*v*k = expected number of detections, where k is number of snapshots (data), v is effective search area, and D is GROUP (not animal) density
        logDgrp[z] <- logDgrp.hyper + logDgrp.err[z]    # expected group density on log scale (uncorrected for g0)
        Dgrp[z] <- exp(logDgrp[z])                      # expected group density on real scale (uncorrected)
        logDgrp.err[z] ~ dnorm(0, logDgrp.tau)          # random DASBR effect for group density (uncorrected)
        Dper1000[z] <- Dgrp[z]*s/g0 * 1000              # animal density, g0 corrected, and rescaled to 'per 1000km2' (note: g0 is the availability corrxn, see below)
      } # Z


      ##### PRIORS FOR DENSITY MODEL #####
      logDgrp.hyper ~ dnorm(0,0.001)              # hyper-parameter for group density, log scale
      logDgrp.tau <- 1/(logDgrp.sig*logDgrp.sig)
      logDgrp.sig ~ dunif(0,4)                    # random effect variance for group density, log scale


      ##### SOME MODEL CHECKS #####
      mean.of.Errs <- mean(logDgrp.err[])         # mean of the estimated random effects (on log scale), should be close to 0
      SD.of.Errs <- sd(logDgrp.err[])             # sd of the estimated random effects, should be similar to logDgrp.sig
      mean.of.logDgrp <- mean(logDgrp[])          # mean of the log-densities across DASBRs, should be similar to logDgrp.hyper


      ##### DERIVED DENSITY ESTIMATE FOR STUDY AREA #####
      Dgrp.mean <- sum(Dgrp[])/nSites  # uncorrected mean density of *groups* in study area (arithmetic mean across DASBR-specific densities)
      D.mean <- Dgrp.mean * s / g0      # animal density = uncorrected group density * mean group size * availability correction
      D.mean.per1000 <- D.mean*1000     # Animal density per 1000 km2
      N <- D.mean * A                   # Population size (A is study area size)
      # Note: In this code, g0 is the availability correction, although the term g0 is used differently in our paper (we call the availability corrxn something else in the paper)


      ##### AN ALTERNATE ESTIMATOR FOR DENSITY - DO NOT USE #####
      Dgrp.mean.II <- exp(logDgrp.hyper + logDgrp.sig*logDgrp.sig/2)                          # group density, uncorrected
      D.mean.per1000.II <- exp(logDgrp.hyper + logDgrp.sig*logDgrp.sig/2) * s / g0 * 1000     # animal density, per 1000km2
      N.II <- D.mean.per1000.II * A/1000                                                      # population size
      # Note: DO NOT USE THIS OUTPUT
      # This estimator seems theoretrically reasonable to me, but results are clearly biased (too high) (due to extreme skew of lognormal?, or due to correlation in the paramemters?)
      # The estimator is the simple formula for finding the mean of a lognormal distribution given lognormal parameters, i.e., mean = exp(mu + sig2/2)
      # A future question is to understand why Dgrp.mean and Dgrp.mean.II don't agree with each other


      ##### ANCILLARY INPUTS TO DENSITY MODELS ABOVE #####
      s ~ dnorm(gsizeMean, tau.s)           # group size.  Mean s comes from Barlow 2016 report
        tau.s <- 1/(gsizeSd * gsizeSd)
        # sd.s <- gsizeSd

      # from EH analysis... number of minutes that a beaker is available for detection...
      timeAvail ~ dnorm(availMean, tau.timeAvail)
        tau.timeAvail <- 1/(sd.timeAvail * sd.timeAvail)
        sd.timeAvail <- availSd
      # Note, the 17.23 = estimate of 15.23 mean avail time + 2 min (for snapshot duration).
      # See Jay's paper on snapshot length for explanation of adding the duration (Barlow 2021, in press)

      # effective g0, the proportion of time between start of successive dives that an animal is available to detection (clicking and oriented properly to hydrophone)
      g0 <- timeAvail / deepDivePeriod
      deepDivePeriod ~ dnorm(ddpMean,ddp.tau)
      ddp.tau <- 1/(ddpSd * ddpSd)
      # ddp.sd <- ddpSd
      # deep dive period (min) is from email from Jay, dated July 23, 2020

    }	# end model

", fill=TRUE, file=file)
    file
}