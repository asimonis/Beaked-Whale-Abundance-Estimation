# robojay functions
library(lubridate)
library(PAMpal)
library(dplyr)
library(mgcv)
library(PAMmisc)
library(readxl)
# Updated 4-20-2022 7pm
formatClickDf <- function(x, stationPattern='(.*)', pascal=FALSE) {
    clicks <- getClickData(x)
    clicks$db <- basename(clicks$db)
    if(pascal) {
        stationPattern <- 'Station-([0-9]{1,2})[-_].*'
    }
    clicks <- select(clicks, all_of(c('UID', 'angle', 'UTC', 'BinaryFile',
                                      'species', 'eventId', 'Latitude', 'Longitude', 'db'))) %>%
        distinct() %>%
        mutate(angle = angle * 180 / pi,
               station = gsub(stationPattern, '\\1', db)) %>%
        group_by(eventId) %>%
        mutate(nClicks = n()) %>%
        ungroup()
    if(pascal) {
        clicks$station <- as.numeric(clicks$station)
        ##### Fix station times
        add7Hours <- clicks$station %in% c(2, 3, 4, 5, 6, 9, 12, 15, 16, 26, 27, 28, 29, 30)
        clicks$UTC[add7Hours] <- clicks$UTC[add7Hours] + 7 * 3600
        add108Hours <- clicks$station == 1
        clicks$UTC[add108Hours] <- clicks$UTC[add108Hours] + 10.8 * 3600
        ######

        clicks$compTime <- floor_date(clicks$UTC, unit='2minute')
        isCont <- clicks$station >= 23
        clicks$compTime[!isCont] <- as.POSIXct(
            gsub('.*_([0-9]{8}_[0-9]{6})\\.pgdf$', '\\1', clicks$BinaryFile[!isCont]),
            format='%Y%m%d_%H%M%S', tz='UTC')

        clicks$compTime[!isCont & add7Hours] <- clicks$compTime[!isCont & add7Hours] + 7 * 3600
        clicks$compTime[!isCont & add108Hours] <- clicks$compTime[!isCont & add108Hours] + 10.8 * 3600
        clicks$diffTime <- as.numeric(difftime(clicks$UTC, clicks$compTime, units='secs'))
        # Jay's done this way, but if have identical binary file name (it happens) is off counts
        # Need to group by station and comp time
        # Bad ex. bin Click_Detector_Click_Detector_Clicks_20160906_195600
        clicks <- clicks %>% group_by(compTime) %>%
            mutate(nCountLT60sec = sum(diffTime < 60),
                   nCountGT60sec = sum(diffTime >= 60)) %>%
            ungroup()
        clicks[c('compTime', 'diffTime')] <- NULL
        clicks <- bind_rows(lapply(split(clicks, clicks$eventId), function(x) {
            x$comment <- ancillary(data[[x$eventId[1]]])$eventComment
            x
        }))
        clicks$eventId <- unname(sapply(clicks$eventId, fmtEventId))
    }
    clicks <- clicks %>% group_by(eventId) %>% mutate(nClicks = n()) %>% ungroup()
    clicks
}

# pilot study only
# jayRenamer <- function(x) {
#     renamer <- c('Angle'='angle', 'eventType' = 'species', 'EventId' = 'eventId',
#                  'lat'='Latitude', 'long'='Longitude', 'ClickDateTime'='UTC',
#                  'Station'='station', 'Amplitude'='dBPP', 'AngleRefractCorrected'='angleCorrected')
#     x <- rename(x, any_of(renamer))
#     x[c("BinaryFile","UID","Angle","ClickDateTime","EventId","Amplitude",
#         "eventType","comment","lat","long","Station","nClicks","nCountLT60sec","nCountGT60sec")]
# }

splitSSPGroups <- function(x, splitBy=NULL, time=Inf, pascal=FALSE) {
    if(pascal) {
        x$UTC <- x[['ClickDateTime']]
    }
    coords <- c('UTC')
    if(!all(coords %in% colnames(x))) {
        warning("Need UTC columns")
        return(x)
    }
    if(is.null(splitBy) ||
       !splitBy %in% colnames(x)) {
        x <- list(x)
    } else {
        x <- split(x, x[[splitBy]])
    }
    tGroup <- 0
    for(i in seq_along(x)) {
        times <- as.numeric(x[[i]]$UTC)
        if(is.infinite(time)) {
            timeCuts <- rep(1, length(times))
        } else {
            timeCuts <- as.numeric(cut(times, seq(from=min(times)-1, to=max(times)+time, by=time)))
        }
        x[[i]]$sspGroup <- tGroup + timeCuts
        tGroup <- max(timeCuts) + tGroup
    }
    x <- bind_rows(x)
    x$sspGroup <- as.character(x$sspGroup)
    if(pascal) {
        x$UTC <- NULL
    }
    x
}

# pascal only
fmtEventId <- function(x) {
    id <- gsub('.*\\.OE([0-9]*)$', '\\1', x)
    if(grepl('Card-', x)) {
        card <- gsub('.*(Card-[A-z]{1,2})_.*', '\\1', x)
        id <- paste0(id, '_', card)
    }
    id
}

createSSPList <- function(x, time=3600*24*1, file=NULL, pascal=FALSE, timeout=240) {
    if(pascal) {
        x <- rename(x, c('Latitude'='lat', 'Longitude'='long', 'UTC'='ClickDateTime'))
    }
    sspCoords <- bind_rows(lapply(split(x, x$sspGroup), function(s) {
        s[which.min(s$UTC)[1], c('UTC','Latitude', 'Longitude', 'sspGroup')]
    }))
    sspList <- createSSP(sspCoords, dropNA=TRUE, timeout=timeout)
    names(sspList) <- sspCoords$sspGroup
    if(!is.null(file)) {
        saveRDS(sspList, file)
    }
    sspList
}

doAngleCorrection <- function(x, sspList, hpDepth=110, animalDepth=1000) {
    ix <- 0
    cat('Calculating raytrace profiles...\n')
    pb <- txtProgressBar(min=0, max=length(sspList), style=3)
    rtAngle <- lapply(sspList, function(s) {
        gamout <- gam(formula=speed ~ s(depth, k=min(30, length(s$speed))), data = s)
        newdata <- data.frame(depth = 0:max(s$depth))
        predicted <- as.numeric(suppressWarnings(predict.gam(gamout, newdata=newdata)))
        newdata$soundSpeed <- predicted
        angles <- 1:87
        # hpDepth <- 110
        # animalDepth <- 1000 + 1
        rt <- raytrace(0, hpDepth, angles, 7, zz=newdata$depth, cc=newdata$soundSpeed, plot = FALSE, progress=FALSE)
        angleOut <- rep(0, length(angles))
        for(i in seq_along(angleOut)) {
            angleOut[i] <- atan((rt$z[[i]][animalDepth+1] - hpDepth)/
                                    rt$x[[i]][animalDepth+1]) * 180 / pi
        }
        ix <<- ix + 1
        setTxtProgressBar(pb, value=ix)
        list(angle=angles, correction=angleOut - angles)
    })


    if(!all(unique(x$sspGroup) %in% names(rtAngle))) {
        warning('Not all groups have corresponding angle corrections, cannot proceed.')
        return(NULL)
    }
    if('Angle' %in% names(x)) {
        x <- rename(x, angle = Angle)
    }
    bind_rows(lapply(split(x, x$sspGroup), function(d) {
        thisAngle <- rtAngle[[x$sspGroup[1]]]
        goodIx <- seq_along(thisAngle$angle)
        tmpAngle <- d$angle - 90
        angleIndex <- round(tmpAngle + 1 - thisAngle$angle[1])
        correction <- rep(0, length(angleIndex))
        for(i in seq_along(correction)) {
            if(angleIndex[i] %in% goodIx) {
                correction[i] <- thisAngle$correction[angleIndex[i]]
            }
        }
        d$correction <- correction
        d$angleCorrected <- tmpAngle + correction
        d$angleCorrected <- ifelse((90 - d$angleCorrected) < 0, 0, 90 - d$angleCorrected)
        d
    }))
}

markSnapshotId <- function(x, snapshot=2, pascal=FALSE) {
    ix <- 0
    x$snapIx <- NULL
    bind_rows(lapply(split(x, x$station), function(s) {
        bind_rows(lapply(split(s, s$eventId), function(e) {
            # isCont <- grepl('/0$', e$dutyCycle[1])
            isCont <- FALSE # testing station 6 issues. Is there a good reason to do other
            # split? I think we can double count either way...
            if(isCont) {
                e$compTime <- floor_date(min(e$UTC), unit=paste0(snapshot, 'min'))
            } else {
                e$compTime <- as.POSIXct(
                    gsub('.*_([0-9]{8}_[0-9]{6})\\.pgdf$', '\\1', e$BinaryFile),
                    format='%Y%m%d_%H%M%S', tz='UTC')
                e$compTime <- min(e$compTime)

                if(pascal) {
                    add7Hours <- e$station %in% c(2, 3, 4, 5, 6, 9, 12, 15, 16, 26, 27, 28, 29, 30)
                    # clicks$UTC[add7Hours] <- clicks$UTC[add7Hours] + 7 * 3600
                    add108Hours <- e$station == 1
                    # clicks$UTC[add108Hours] <- clicks$UTC[add108Hours] + 10.8 * 3600
                    ######

                    # clicks$compTime <- floor_date(clicks$UTC, unit='2minute')
                    # isCont <- clicks$station >= 23
                    # clicks$compTime[!isCont] <- as.POSIXct(
                    #     gsub('.*_([0-9]{8}_[0-9]{6})\\.pgdf$', '\\1', clicks$BinaryFile[!isCont]),
                    #     format='%Y%m%d_%H%M%S', tz='UTC')
                    #
                    e$compTime[add7Hours] <- e$compTime[add7Hours] + 7 * 3600
                    e$compTime[add108Hours] <- e$compTime[add108Hours] + 10.8 * 3600
                }
            }

            e$diffMin <- as.numeric(difftime(e$UTC, e$compTime, units='mins'))
            e$snapIx <- floor(e$diffMin / snapshot)
            e$snapTime <- e$compTime + e$snapIx * 60 * snapshot
            e$snapIx <- e$snapIx + 1 + ix
            ix <<- max(e$snapIx)
            # e$compTime <- NULL
            # e$diffMin <- NULL
            e
        }))
    }))
}

formatEventSummary <- function(x, recorderInfo=NULL, snapshot=2, pascal=FALSE) {
    evCol <- c('EventId', 'eventId')[c('EventId', 'eventId') %in% colnames(x)]
    latCol <- c('lat', 'Latitude')[c('lat', 'Latitude') %in% colnames(x)][1]
    longCol <- c('long', 'Longitude')[c('long', 'Longitude') %in% colnames(x)][1]
    # this is currently split by binary to do averaging because each binary is 2 minutes
    # which is our snapshot.
    result <- addRecorderInfo(x, info=recorderInfo, pascal=pascal)
    result <- markSnapshotId(result, snapshot=snapshot, pascal=pascal)
    # mark snapshot ids.
    # result <- bind_rows(lapply(split(x, x$BinaryFile), function(b) {
    # bind_rows(lapply(split(b, b[[evCol]]), function(e) {
    result <- bind_rows(lapply(split(result, result$snapIx), function(e) {
        # e$fileDateTime = as.POSIXct(
        #     gsub('.*_([0-9]{8}_[0-9]{6})\\.pgdf$', '\\1', e$BinaryFile),
        #     format='%Y%m%d_%H%M%S', tz='UTC')
        e$fileDateTime <- e$snapTime[1]
        e[[latCol]] <- e[[latCol]][1]
        e[[longCol]] <- e[[longCol]][1]
        if('UTC' %in% colnames(e)) {
            e$UTC <- min(e$UTC)
        }
        angles <- e$angleCorrected
        histOut <- hist(angles,breaks=seq(0,181,1),plot=FALSE)
        findMaxFreq <- which.max(histOut$counts[c(1:88,92:181)])      #find modal value, exlc values around 90 deg.
        if(findMaxFreq > 88) {
            findMaxFreq <- findMaxFreq + 3
        }
        modalAngle <- histOut$mids[findMaxFreq]
        #Calculate mean detection angle for angles less than 88 deg if Modal Angles is less than 88 det
        meanAngle= NA
        useAngles <- angles < 88 & e$correction != 0
        if(modalAngle < 88 & any(useAngles)) {
            meanAngle= atan(mean(tan(angles[useAngles]*pi/180),na.rm=TRUE)) * 180 / pi  #average angles is within 2 deg of mode
        }
        # if(is.na(meanAngle)) {
        #     plot(angles)
        #     browser()
        # }
        # if(!is.na(meanAngle) && abs(modalAngle-meanAngle) > 3) {
        #     plot(x=e$angleCorrected)
        # }
        e$modalAngle <- modalAngle
        e$meanAngle <- meanAngle
        e$nClicks <- nrow(e)

        # select(e, any_of(c('UTC', 'fileDateTime', 'BinaryFile', evCol, 'species','eventType',
        #                    'meanAngle', 'modalAngle', 'station', 'nClicks', latCol, longCol,
        #                    'comment', 'nCountLT60sec', 'nCountGT60sec', 'db', 'snapIx', 'compTime')))
        select(e, any_of(c('UTC', 'fileDateTime', evCol, 'species','eventType',
                           'meanAngle', 'modalAngle', 'station', 'nClicks', latCol, longCol,
                           'comment', 'db', 'snapIx', 'dutyCycle', 'recorder')))

        # }))
    }))
    result <- arrange(result, station, fileDateTime)
    result <- distinct(result)
    # browser()
    if(pascal) {
        # add7Hours <- result$station %in% c(2, 3, 4, 5, 6, 9, 12, 15, 16, 26, 27, 28, 29, 30)
        # result$fileDateTime[add7Hours] <- result$fileDateTime[add7Hours] + 7 * 3600
        # add108Hours <- result$station == 1
        # result$fileDateTime[add108Hours] <- result$fileDateTime[add108Hours] + 10.8 * 3600
        result$recorder <- 'ST4300'
        result$recorder[result$station %in% c(7, 10, 13, 17)] <- 'SM3M'
        result$recorder[result$station %in% c(8, 11, 14, 18)] <- 'SM2Bat'
        # result$MultipleAngles <- grepl('MA', result$comment)
        # result$Analyst= substr(result$comment,1,1)
        # splitComment <- strsplit(result$comment, ',')
        # ici <- sapply(splitComment, function(s) s[7])
        # result$ICI <- as.numeric(ici)
    }
    result$deltaMin <- NA
    result$deltaMin[2:nrow(result)] <- as.numeric(difftime(
        result$fileDateTime[2:nrow(result)],
        result$fileDateTime[1:(nrow(result)-1)],
        units='mins'))
    result$deltaMin[result$deltaMin < 0] <- NA
    result
}

lognorm<- function(n, mean, sd) {
    meanlog <- log(mean^2 / sqrt(sd^2 + mean^2))
    sdlog<- sqrt(log(1 + (sd^2 / mean^2)))
    rlnorm(n, meanlog, sdlog)
}

createSimLikeFun <- function(nSim=1e6, model=c('HN', 'C_HN', 'HR'), depthDistr=c('log-normal'),
                             meanDepth=1217,sdDepth=354, sdAngleErr=0.59,
                             mainTitle='', truncAngle=NULL, hpDepth=110, truncDist=4e3,
                             par1Lim=c(1e3, 3e3), par2Lim=c(-2, 28)) {
    # maxHRange<- 8000
    maxHRange <- truncDist * 2
    if(is.null(truncAngle)) {
        truncAngle<- atan(truncDist/(meanDepth-hpDepth))*180/pi
    }
    hRange <- sqrt((maxHRange^2)*runif(nSim))
    depth <- if(depthDistr == "log-normal") {
        lognorm(nSim,mean=(meanDepth-hpDepth),sd=sdDepth)   # distribution depth (m)
    } else if (depthDistr == "normal") {
        rnorm(nSim,mean=(meanDepth-hpDepth),sd=sdDepth)     # distribution depth (m)
    } else {
        cat("  ERROR, need to specify the depth distribution function as normal or log-normal")
    }
    depth[depth < 300] <- 300

    sRange <- sqrt(hRange^2 + depth^2)
    angle <- abs(rnorm(nSim, mean = atan(hRange/depth)*180/pi, sd=sdAngleErr))
    doTrunc <- angle < truncAngle
    sRange <- sRange[doTrunc]
    angle <- angle[doTrunc]
    # dt <- data.table(iAngle = c(ceiling(Angle), 1:91))
    dt <- data.frame(iAngle = factor(c(ceiling(angle), 1:91), levels=1:91))
    # setkey(dt, iAngle)
    # negRangeSq <- -(SRange^2)/2
    model <- match.arg(model)
    function(detAngle) {
        detAngle <- detAngle[detAngle < truncAngle & !is.na(detAngle)]
        detAngle <- ceiling(detAngle)
        switch(model,
               'C_HN' = {
                   probFun <- function(param) {
                       if(param[1] < par1Lim[1] ||
                          param[1] > par1Lim[2] ||
                          param[2] < par2Lim[1] ||
                          param[2] > par2Lim[2]) {
                           return(rep(exp(-1000), length(sRange)))
                       }
                       1- (1-exp(-.5*(sRange/param[1])^2))^(exp(param[2]))
                   }
               },
               'HN' = {
                   probFun <- function(param) {
                       if(param[1] < par1Lim[1] ||
                          param[1] > par1Lim[2]) {
                           return(rep(exp(-1000), length(sRange)))
                       }
                       exp(-.5*(sRange/param[1])^2)
                   }
               },
               'HR' = {
                   probFun <- function(param) {
                       if(param[1] < par1Lim[1] ||
                          param[1] > par1Lim[2] ||
                          param[2] < par2Lim[1] ||
                          param[2] > par2Lim[2]) {
                           return(rep(exp(-1000), length(sRange)))
                       }
                       1 - exp(-(sRange/param[1])^(-param[2]))
                   }
               }
        )
        function(param, like=TRUE) {
            probDet <- probFun(param)
            dt$probDet <- c(probDet, rep(0, 91))
            sumprob <- group_by(dt, iAngle) %>% summarise(sum=sum(probDet)) %>% .$sum
            # dt[, prob := probDet]
            # sumprob <- dt[, .(sum(prob)), keyby=iAngle][[2]]
            # dt[, prob:= NULL]
            sumprob <- sumprob / sum(sumprob)
            # browser()
            if(like) {
                return(sum(log(sumprob[detAngle])))
            }
            hRange <- hRange[doTrunc]
            depth <- depth[doTrunc]
            # calculate average prob of detection in 100 categories of horizontal range
            probDetHRange <- rep(NA,100)
            dx <- maxHRange/100
            HRindex <- ceiling(100*hRange/maxHRange)
            HR <- (1:100 - 0.5)*dx
            HRindex <- factor(HRindex, levels=1:100)
            probDetHRange <- sapply(split(probDet, HRindex), function(p) mean(p, na.rm=TRUE))
            
            integralRangeTimesDetProb <- sum(HR[HR<=truncDist] * probDetHRange[HR<=truncDist],na.rm = TRUE)*dx
            EDR <- sqrt(2*integralRangeTimesDetProb)
            # calculate the EDR for 39 independent 2-minute detection opportunities for dive
            probDet39X <- 1 - (1-probDetHRange)^39
            integralRangeTimesDetProb <- sum(((1:100 - 0.5)*dx) * probDet39X,na.rm = TRUE)*dx
            EDR_39X <- sqrt(2*integralRangeTimesDetProb)
            
            list(angleDensity=sumprob,probDetHRange=probDetHRange,EDR=EDR,angle=angle,
                 probDet=probDet,EDR_39X=EDR_39X, sRange=sRange)
        }
    }
}

newEstDetFunction <- function(eventAngles,
                              species = NULL,
                              subsetCol = 'recorder',
                              meanDepth = 1217,
                              sdDepth = 354,
                              hpDepth = 110,
                              truncDist = 4000,
                              sdAngleErr = 0.59,
                              depthDistr = 'log-normal',
                              range = c(750, 3000),
                              model = c('C_HN', 'HN', 'HR'),
                              outDir= '.',
                              doJackknife = FALSE,
                              jk_nSamples = 20,
                              nSim = 1e7,
                              verbose=FALSE,
                              progress=TRUE,
                              simSpread=FALSE,
                              ...) {

    truncAngle<- atan(truncDist/(meanDepth-hpDepth))*180/pi

    eventAngles<- eventAngles[(eventAngles$nClicks>2),]
    model <- match.arg(model)
    # Limit sample to SpeciesID
    if(!is.null(species)) {
        spCol <- c('species', 'eventType')[c('species', 'eventType') %in% colnames(eventAngles)][1]
        eventAngles<- eventAngles[eventAngles[[spCol]] == species,]
    }
    subsets <- 'All'
    if(!is.null(subsetCol)) {
        subsets <- c(subsets, unique(eventAngles[[subsetCol]]))
    }
    outs <- vector('list', length=length(subsets))
    names(outs) <- subsets
    if(progress) {
        cat('Running model on subsets...\n')
        nRuns <- ifelse(length(subsets) == 2, 1, length(subsets))
        pb <- txtProgressBar(min=0, max=nRuns*(1 + doJackknife * jk_nSamples), style=3)
        pbIx <- 0
    }
    for(subsetName in subsets) {
        mainTitle<- paste("PASCAL",subsetName,model,sep=" ")
        # Send output to a file
        outputFile<- paste(model,"DetFunctFit","n",nSim,truncDist,subsetName,".txt",sep="_")
        outputFile <- file.path(outDir, outputFile)
        if(subsetName == 'All') {
            obsDetAngles <- eventAngles$meanAngle
        } else {
            obsDetAngles <- eventAngles$meanAngle[eventAngles[[subsetCol]] == subsetName]
        }

        # Observed detection angles are valid average angles
        obsDetAngles<- obsDetAngles[!is.na(obsDetAngles)]
        if(simSpread) {
            simAngles <- vector('list', length=length(obsDetAngles))
            for(s in seq_along(simAngles)) {
                if(runif(1) > .33) {
                    simAngles[[s]] <- obsDetAngles[s]
                    next
                }
                simAngles[[s]] <- rnorm(3, mean=obsDetAngles[s], sd=5)
            }
            obsDetAngles <- unlist(simAngles)
        }
        
        # plot distribution of detection angles before truncation
        if(verbose) {
            cat(" Tot number of files w/ ",species," and 3 or more clicks = ",length(obsDetAngles))
            # if (SubSetName == "All Recorders") { upperFreq<-80 }else{upperFreq<-40}
            histout<-hist(obsDetAngles,breaks=50,col="black",xlab="Detection Angle (deg)",
                          xlim=c(0,100),main=mainTitle)
            lines(c(truncAngle,truncAngle),c(0,80),lty="dashed",col="white",lwd=2)
            maxcount<- histout$breaks[which.max(histout$counts)]
            hDistAtMaxCount<- (meanDepth-hpDepth) * tan(maxcount*pi/180)
            cat(" Maximum count at =",maxcount,"deg, corresponding to a Hdist of ",hDistAtMaxCount,"m \n")
            # truncate observed detection angles
            pctTrunc<- 100*sum(obsDetAngles>truncAngle,na.rm=T)/sum(!is.na(obsDetAngles))
            cat(" Percentage of detection angles truncated =",pctTrunc,"% \n")
        }
        obsDetAngles<- obsDetAngles[obsDetAngles < truncAngle]
        nSample<- length(obsDetAngles)
        par2Lim <- switch(model,
                     'C_HN' = c(-2, 28),
                     'HN' = 0,
                     'HR' = c(.01, 50)
        )
 
        thisLikeBase <- createSimLikeFun(nSim=nSim, model=model, depthDistr=depthDistr,
                                         meanDepth=meanDepth, sdDepth=sdDepth, sdAngleErr=sdAngleErr,
                                         mainTitle=mainTitle, truncAngle=truncAngle, hpDepth=hpDepth,
                                         truncDist=truncDist, par1Lim=range(range), par2Lim=par2Lim)
        thisLikeFun <- thisLikeBase(obsDetAngles)
        opResult <- switch(
            model,
            'C_HN' = optim(c(1e3, 1), fn=thisLikeFun, method='Nelder-Mead', hessian=FALSE,
                           control=list(fnscale=-1, parscale=c(mean(range(range)), mean(range(par2Lim))))),
            'HN' =  optim(c(1e3), fn=thisLikeFun, method='Brent', upper=truncDist, lower=500, hessian=FALSE,
                          control=list(fnscale=-1)),
            'HR' = optim(c(1e3, 1), fn=thisLikeFun, method='Nelder-Mead', hessian=FALSE,
                         control=list(fnscale=-1, parscale=c(mean(range(range)), mean(range(par2Lim)))))
        )
        
        outVals <- thisLikeFun(opResult$par, like=FALSE)

        maxL_EDR<- outVals$EDR
        maxL_EDR39X<- outVals$EDR_39X

        hRangeDetProbDF<- data.frame(hRange=(1:100)*80,detProb=outVals$probDetHRange)
        # plot cumulative distribution of angles from observations and simulation

        maxL_ESA<- pi*maxL_EDR^2  #effective survey area
        inclAngle<- outVals$angle[runif(length(outVals$probDet))<outVals$probDet]
        KS_out<- ks.test(obsDetAngles,inclAngle,exact=FALSE)  #test of observed and predicted angles
        if(verbose) {
            cdf_Angles<- sapply(1:91,function(x) mean(inclAngle<=x))
            cdf_obsDetAngless<- sapply(1:91,function(x) mean(obsDetAngles<=x))
            plot(1:91,cdf_Angles,type="l",lwd=3,main=mainTitle,xlab="Detection Angle (deg.)",
                 ylab="Cumulative Distribution")
            lines(1:91,cdf_obsDetAngless,col="red",lwd=3)


            # plot((1:100)*80,outVals$probDetHRange,type="l",lwd=5,main=mainTitle,xlab="Horizontal Range (m)",
            #      ylab="Probability of Detection",ylim=c(0,1.05),xlim=c(0,6000))
            # plot((1:100)*80/1000,outVals$probDetHRange,type="l",lwd=5,main=NULL,xlab="Horizontal Range (km)",
            #      ylab="Probability of Detection",ylim=c(0,1.05),xlim=c(0,5))
            plotOneDetFun(opResult$par, model=model, truncDist=truncDist)
            cat(" Sample size of detection angles (#snapshots):",nSample,"\n")
            cat(" Detection function for: ",mainTitle,"\n")
            cat(" Maximum likelihood estimate of EDR:",maxL_EDR/1000,"km \n")
            cat(" Maximum likelihood estimate of ESA:",maxL_ESA/(1000^2),"km2 \n")
            cat(" ML goodness of fit from KS stat, p=",KS_out$p.value,"\n")
            cat(" The log-likelihood of ML estimates is",opResult$value,"\n")
            cat(" AIC =",-2*(opResult$value-length(opResult$par)),"\n")
            cat(" Prob detection at zero horiz range (gzero)=",hRangeDetProbDF$detProb[1])
        }
        if(progress) {
            pbIx <- pbIx + 1
            setTxtProgressBar(pb, value=pbIx)
        }
        thisOut <- list(maxLikeParam=opResult$par,nSample=nSample, mainTitle=mainTitle, maxL_EDR=maxL_EDR,
                        maxL_ESA=maxL_ESA, KS_prob=KS_out$p.value, MaxLik=opResult$value,
                        AIC = -2*(opResult$value-length(opResult$par)),
                        gzero = hRangeDetProbDF$DetProb[1],
                        hRangeDetProbDF=hRangeDetProbDF,
                        angle=obsDetAngles,
                        likeFun=thisLikeBase,
                        model=model)
        if(doJackknife) {
            jk_EDR<- array(1)
            jk_EDR_39X<- array(1)
            nAngles<- length(obsDetAngles)
            jk_sample<- ceiling((1:nAngles) * jk_nSamples / nAngles)

            for (j in 1:jk_nSamples) {
                gc()
                jkStart <- Sys.time()
                jk_obsDetAngles<- obsDetAngles[!(jk_sample == j)]
                jkLikeFun <- thisLikeBase(jk_obsDetAngles)
                ##########################################################################
                # get likelihood surface data for observed angles by simulated likelihood
                ##########################################################################
                jkOptim <- if(model == 'C_HN') {
                    optim(c(1e3, 1), fn=jkLikeFun, method='Nelder-Mead', hessian=FALSE,
                          control=list(fnscale=-1, parscale=c(mean(range(range)), mean(range(par2Lim)))))
                } else {
                    optim(c(1e3, 1), fn=jkLikeFun, method='Brent', upper=4e3, lower=500, hessian=FALSE, control=list(fnscale=-1))
                }
                # find maximum likelihood value from discrete values of detection function parameters
                jkOuts <- jkLikeFun(jkOptim$par, like=FALSE)
                if(verbose) {
                    cat('\nJK optim #', j, 'count:', jkOptim$counts[1],
                        'Time:', round(as.numeric(difftime(Sys.time(), jkStart)), 1))
                    cat(" Jackknife sample:",j,"\n")
                    cat(" Sample size of detection angles (#snapshots):",length(jk_obsDetAngles),"\n")
                    cat(" Maximum log likelihood from likelihood surface data is:",jkOptim$value,
                        " with Param1=",jkOptim$par[1],"and Param2=",jkOptim$par[2],"\n")
                    cat(" Maximum likelihood estimate of EDR:",jkOuts$EDR,"\n")
                }

                jk_EDR[j]<- jkOuts$EDR
                jk_EDR_39X[j]<- jkOuts$EDR_39X
                if(progress) {
                    pbIx <- pbIx + 1
                    setTxtProgressBar(pb, value=pbIx)
                }
            }
            ################################################################################
            # Summary output for snapshot (2-min) detection probability
            nJK <- length(jk_EDR)
            mean_jk_EDR <- mean(jk_EDR)
            thisOut$jk_EDR <- jk_EDR
            thisOut$mean_jk_EDR <- mean_jk_EDR
            JackKnife_var <- ((nJK-1)/nJK) * sum((jk_EDR-mean_jk_EDR)^2)
            thisOut$jackKnife_var <- JackKnife_var
            thisOut$jackKnife_cv <- sqrt(JackKnife_var)/mean_jk_EDR

            JK_ESA<- pi*jk_EDR^2  #effective survey area
            mean_JK_ESA<- mean(JK_ESA)
            thisOut$mean_jk_ESA <- mean_JK_ESA
            Jackknife_varESA<- ((nJK-1)/nJK) * sum((JK_ESA-mean_JK_ESA)^2)
            thisOut$jackknife_varESA <- Jackknife_varESA
            thisOut$jackKnife_cvESA <- sqrt(Jackknife_varESA)/mean_JK_ESA
            if(verbose) {
                cat(" Maximum likelihood EDR (1X)=",maxL_EDR,"m \n")
                cat(" Jackknife mean EDR=",mean_jk_EDR,"  se =",sqrt(JackKnife_var),"  cv =",
                    sqrt(JackKnife_var)/mean_jk_EDR,"\n")
                cat(" Jackknife mean ESA=",mean_JK_ESA,"  se =",sqrt(Jackknife_varESA),"  cv =",
                    sqrt(Jackknife_varESA)/mean_JK_ESA, "\n")
            }
            ################################################################################
            # Summary output for dive detection probability (estimated a 39 independent snapshots)

            nJK<- length(jk_EDR_39X)
            mean_JK_EDR_39X<- mean(jk_EDR_39X)
            thisOut$jk_EDR_39x <- jk_EDR_39X
            thisOut$mean_jk_EDR_39X <- mean_JK_EDR_39X
            JackKnife_var<- ((nJK-1)/nJK) * sum((jk_EDR_39X-mean_JK_EDR_39X)^2)
            thisOut$jackKnife_var_39x <- JackKnife_var
            thisOut$jackKnife_cv_39x <- sqrt(JackKnife_var)/mean_JK_EDR_39X

            JK_ESA_39X<- pi*jk_EDR_39X^2  #effective survey area
            mean_JK_ESA_39X<- mean(JK_ESA_39X)
            thisOut$mean_jk_ESA_39X <- mean_JK_ESA_39X
            Jackknife_varESA_39X<- ((nJK-1)/nJK) * sum((JK_ESA_39X-mean_JK_ESA_39X)^2)
            thisOut$jackknife_varESA_39X <- Jackknife_varESA_39X
            thisOut$jackKnife_cvESA_39x <- sqrt(Jackknife_varESA_39X)/mean_JK_ESA_39X
            if(verbose) {
                cat(" Maximum likelihood dive EDR (39X)=",maxL_EDR39X,"m \n")
                cat(" Jackknife mean dive EDR (39X)=",mean_JK_EDR_39X,"  se =",sqrt(JackKnife_var),"  cv =",
                    sqrt(JackKnife_var)/mean_JK_EDR_39X,"\n")
                cat(" Jackknife mean dive ESA (39X)=",mean_JK_ESA_39X,"  se =",sqrt(Jackknife_varESA_39X),"  cv =",
                    sqrt(Jackknife_varESA_39X)/mean_JK_ESA_39X, "\n")
            }
        }
        outs[[subsetName]] <- thisOut
        # dont repeat calcs if we only have 1 subset since we always have "All"
        if(length(subsets) == 2) {
            outs[[2]] <- thisOut
            break
        }
    } # subset loop end
    outs
}

createDetectionHistory <- function(ea, maxTime=61.999, meanDepth =1217, sdDepth=354,
                                   hpDepth = 110, truncDist = 4000, plot=FALSE, snapshot=2) {
    truncAngle <- atan(truncDist/(meanDepth-hpDepth))*180/pi  #actual truncation angle based on nominal truncation dist. and mean depth

    if(is.character(ea$fileDateTime)) {
        ea$fileDateTime <- as.POSIXct(ea$fileDateTime, "%Y-%m-%d %H:%M:%S", tz='UTC')
    }
    # Why is this done and then undone in Jeff's code
    ea$meanAngle <- 180 - ea$meanAngle
    ea <- ea[ea$nClicks >= 3,]
    ea <- ea[!is.na(ea$meanAngle),]
    ea$station <- as.character(ea$station)
    # limit sample to Ziphius
    ea$deltaMin[is.na(ea$deltaMin)]= -9999
    nLines <- length(ea$deltaMin)

    # plot distribution of detection angles and truncation angle
    if(plot) {
        op <- par(xaxs="i",yaxs="i",cex.axis=1.3,cex.lab=1.5,font=2,font.lab=2)  #format plots
        on.exit(par(op))
        hist(180-ea$meanAngle,xlab="Detection Angle (deg)",ylim=c(0,160),breaks=seq(0,100,5),main=NULL)
        lines(c(truncAngle,truncAngle),c(0,200),lty="dashed")
    }

    # write ZC event file for input to density analysis
    ##### INPUT FOR JEFFS FUNCTiON encounterRate-DataPrep####
    # write.csv(ea,file="ZcEventFilesWithAngles_drifts_1-22 - CorrectedAngles v2.csv", row.names=FALSE)

    # limit samples to those that will be used for abundance estimation (truncation distance and minimum click count)
    ea_trunc <- ea[ea$meanAngle>=(180-truncAngle),]
    # estimate the number of snapshots with detections for each drift w/in trunc angle
    eventsPerStation <- group_by(ea_trunc, station) %>%
        summarise(numEvents = n())
    # making dcycle stuff
    dcDf <- distinct(ea[c('station', 'dutyCycle')])
    dcMap <- lapply(dcDf$dutyCycle, function(x) {
        as.numeric(strsplit(x, '/')[[1]])
    })
    names(dcMap) <- as.character(dcDf$station)
    minDc <- min(sapply(dcMap, sum))
    nCols <- max(sapply(dcMap, function(x) {
        dcToCols(x, maxTime = maxTime, snapshot=snapshot)
    }))

    colMap <- lapply(dcMap, function(x) {
        dcToColMap(x, snapshot=snapshot, ncol=nCols)
    })
    # Search through Ziphius events and assign a unique number for each dive (using MaxTime limit)
    ea <- arrange(ea, station, fileDateTime)
    ea$diveNum= NA
    ea$diveNum[1]= 1
    ea$diveTime= NA
    ea$diveTime[1]= 0
    diveStart= ea$fileDateTime[1]
    dNum= 1
    for(i in 2:length(ea$diveNum)) {
        # ea$diveTime[i]= as.numeric(difftime(ea$fileDateTime[i],diveStart,units="mins"))
        ea$diveTime[i]= as.numeric(difftime(ea$UTC[i],diveStart,units="mins"))
        # browser()
        if((ea$station[i] != ea$station[i-1]) ||
           (ea$diveTime[i] > maxTime) ||
           (ea$diveTime[i] < 0)) {
            dNum= dNum + 1
            ea$diveTime[i]= 0
            diveStart= ea$fileDateTime[i]
        }
        ea$diveNum[i]= dNum
    }
    cat(" total number of dives:",dNum,"\n")
    # browser()
    # Create Detection History dataframe
    nDives= max(ea$diveNum)
    # is this div by 2 because 2 minutes?? And then is max time 61.999 so that div/2 is lower?
    # detectionHistory= array(0,dim=c(nDives,1+maxTime/snapshot))
    detectionHistory <- array(0, dim=c(nDives, nCols))
    recorder= rep(NA,nDives)
    station= rep(NA,nDives)
    time= rep(NA,nDives)
    eventId= rep(NA,nDives)
    ea$filesPerDive= NA
    deltaTime= rep(NA,nDives)
    startAngle= rep(NA,nDives); endAngle= rep(NA,nDives)
    nAngles= 0; angles=rep(0,2); angleTimes= rep(0,2)
    # create dectection history for each recognized dive
    for (idive in 1:nDives) {
        # browser()
        diveEvents <- ea[ea$diveNum == idive,]
        n <- length(diveEvents$meanAngle)
        recorder[idive] <- as.character(diveEvents$recorder[1])
        station[idive] <- diveEvents$station[1]
        thisDcMap <- colMap[[station[idive]]]
        thisTimes <- cumsum(thisDcMap)
        eventId[idive] <- as.character(diveEvents$eventId[1])
        time[idive] <- as.character(diveEvents$fileDateTime[1])
        deltaTime[idive] <- difftime(diveEvents$fileDateTime[n],diveEvents$fileDateTime[1],units='min')
        startAngle[idive] <- diveEvents$meanAngle[1]
        endAngle[idive] <- diveEvents$meanAngle[n]
        for (i in 1:n) {
            elapsedMin= diveEvents$diveTime[i]
            if (elapsedMin < maxTime) {
                # cix <- 1+(elapsedMin/snapshot)
                cix <- min(which(thisTimes > elapsedMin))
                detectionHistory[idive, cix]= diveEvents$meanAngle[i]
                nAngles= nAngles + 1
                angles[nAngles]= diveEvents$meanAngle[i]
                angleTimes[nAngles]= difftime(diveEvents$fileDateTime[i],diveEvents$fileDateTime[1],units='min')
            }
        }
    }
    nEvOff <- 0
    whichDataOff <- numeric(0)
    whichEvOff <- character(0)
    nDetOff <- 0
    debugDf <- list()
    # for instruments on duty cycle, indicate when detections were not possible with -1
    for(i in 1:nrow(detectionHistory)) {
        stationOff <- markOffCycle(ea$dutyCycle[ea$station == station[i]][1], snapshot=snapshot, ncol=nCols, maxTime=maxTime)
        # if(i == 8) browser()
        # duty cycle might be offset if multiple snapshots in one on cycle depending on
        # which one first deteciton is in

        if(any(stationOff)) {
            nOn <- min(which(stationOff)) - 1
            if(nOn > 1) {
                # cat('\n', i)
                # browser()
                # move the duty ycle marks instead of the data
                tryOff <- lapply(0:(nOn-1), function(n) {
                    markOffCycle(ea$dutyCycle[ea$station == station[i]][1], snapshot=snapshot, ncol=nCols, maxTime=maxTime, shift=n)
                })
                nDropped <- sapply(tryOff, function(o) {
                    sum(detectionHistory[i, o] != 0)
                })
                stationOff <- tryOff[[which.min(nDropped)]]

                # tryOn <- lapply(1:nOn, function(n) {
                #     c(rep(0, n-1), detectionHistory[i, 1:(ncol(detectionHistory)-n+1)])
                # })
                # nDropped <- sapply(tryOn, function(t) {
                #     sum(t[stationOff] != 0)
                # })
                # detectionHistory[i, ] <- tryOn[[which.min(nDropped)]]
            }
            # mark how many are getting coverd up as "off cycle"
            hasData <- detectionHistory[i, ] != 0
            thisDataOff <- hasData & stationOff
            
            # thisDataOff <- sum(detectionHistory[i, stationOff] != 0)
            if(any(thisDataOff)) {
                offDf <- ea[ea$diveNum == i, ]
                offDf$markedOff <- FALSE
                offEvIx <- which(which(hasData) %in% which(thisDataOff))
                offEv <- ea$eventId[ea$diveNum == i][offEvIx]
                offDf$markedOff[offEvIx] <- TRUE
                nDetOff <- nDetOff + sum(ea$nClicks[ea$diveNum == i][offEvIx])
                whichDataOff <- c(whichDataOff, i)
                whichEvOff <- c(whichEvOff, unique(offEv))
                debugDf[[length(debugDf)+1]] <- offDf
            }
            nEvOff <- nEvOff + sum(thisDataOff)
            detectionHistory[i, stationOff] <- -1
        }
    }
    if(nEvOff > 0) {
        # warning(nDetOff, ' detections were marked as "off duty cycle"',
                # ' in dive number ', paste0(whichDataOff, collapse=', '))
        warning(nDetOff, ' detections marked as "off duty cycle"',
                ' from ', nEvOff, ' events. Check $debug in the output.')
        
    }
    # create dataframe from detection history matrix and save as csv file
    # detectionHistoryOut= data.frame(station,time,eventId,recorder,detectionHistory)
    detHistDf <- data.frame(station, time, eventId, recorder)
    detectionHistory <- detectionHistory[!is.na(detHistDf$recorder), ]
    detHistDf <- detHistDf[!is.na(detHistDf$recorder), ]
    # detectionHistoryOut= list(dhDf = detHistDf, dhMat = detectionHistory)
    #### OUTPUT FOR JEFF EH-dataPrep ####
    # write.csv(detectionHistoryOut,file="DetectionHistory.csv", row.names=FALSE)
    if(nEvOff > 0) {
        debugDf <- bind_rows(debugDf)
        debugDf <- debugDf[c('UTC', 'eventId', 'nClicks', 'dutyCycle', 'diveNum', 'diveTime', 'markedOff')]
        # debugDf$markedOff <- debugDf$eventId %in% whichEvOff
    } 
        
    if(plot) {
        # determine duration of encounter w/o angle truncation
        detectDuration= rep(0,nDives)
        for (idive in 1:nDives) {
            detectDuration[idive]= (max(which(detectionHistory[idive,]>0)) - min(which(detectionHistory[idive,]>0)) + 1) * 2
        }
        hist(detectDuration[recorder=="SM3M"],breaks=seq(0,80,2),xlab="SM3M Encounter Duration (min)",
             main=paste("Mean =",mean(detectDuration[recorder=="SM3M"]),"  SD=",sd(detectDuration[recorder=="SM3M"])))
        hist(detectDuration[recorder=="SM2Bat"],breaks=seq(0,80,2),xlab="SM2Bat Encounter Duration (min)",
             main=paste("Mean =",mean(detectDuration[recorder=="SM2Bat"]),"  SD=",sd(detectDuration[recorder=="SM2Bat"])))
        hist(detectDuration[recorder=="ST4300"],breaks=seq(0,80,2),xlab="ST4300 Encounter Duration (min)",
             main=paste("Mean =",mean(detectDuration[recorder=="ST4300"]),"  SD=",sd(detectDuration[recorder=="ST4300"])))

        # determine duration of encounter w/ angle truncation
        detectDurationTrunc= rep(0,nDives)
        detectionHistoryTrunc= detectionHistory
        detectionHistoryTrunc[(180-detectionHistoryTrunc) > truncAngle]= 0
        for (idive in 1:nDives) {
            nDet= sum(detectionHistoryTrunc[idive,]>0)
            if (nDet > 0) {
                detectDurationTrunc[idive]= duration= (max(which(detectionHistoryTrunc[idive,]>0)) - min(which(detectionHistoryTrunc[idive,]>0)) + 1) * 2
            }
        }
        hist(detectDurationTrunc[recorder=="SM3M" & detectDurationTrunc>0],breaks=seq(0,80,2),xlab="SM3M Encounter Duration (min)",
             main=paste("Mean =",mean(detectDurationTrunc[recorder=="SM3M" & detectDurationTrunc>0],na.rm=T),"  SD=",
                        sd(detectDurationTrunc[recorder=="SM3M" & detectDurationTrunc>0],na.rm=T)))
        hist(detectDurationTrunc[recorder=="SM2Bat" & detectDurationTrunc>0],breaks=seq(0,80,2),xlab="SM2Bat Encounter Duration (min)",
             main=paste("Mean =",mean(detectDurationTrunc[recorder=="SM2Bat" & detectDurationTrunc>0],na.rm=T),"  SD=",
                        sd(detectDurationTrunc[recorder=="SM2Bat" & detectDurationTrunc>0],na.rm=T)))
        hist(detectDurationTrunc[recorder=="ST4300" & detectDurationTrunc>0],breaks=seq(0,80,2),xlab="ST4300 Encounter Duration (min)",
             main=paste("Mean =",mean(detectDurationTrunc[recorder=="ST4300" & detectDurationTrunc>0],na.rm=T),"  SD=",
                        sd(detectDurationTrunc[recorder=="ST4300" & detectDurationTrunc>0],na.rm=T)))

        plot(deltaTime,(startAngle-endAngle),type="p")
        abline(h=0)
        plot(startAngle,endAngle,type="p",xlim=c(0,90),ylim=c(0,90))
        lines(c(0,90),c(0,90),col="red")
        plot(startAngle,(startAngle-endAngle),type="p")
        abline(h=0)
    }
    # Angles= 180 - Angles
    # plot(Angles[1:(nAngles-1)],Angles[2:nAngles],type="l")
    # plot(AngleTimes,Angles,type="l",xlab="Elapsed Time (min)")
    list(ea=ea,
         dhDf=detHistDf,
         dhMat = detectionHistory,
         dcMap = colMap,
         debug = debugDf)
}

dcToCols <- function(x, maxTime=61.999, snapshot=2) {
    # if(x[1] %% snapshot != 0) {
    #     warning('Recording "on" interval must be multiple of snapshot length')
    #     return(NULL)
    # }
    nDc <- ceiling(maxTime / sum(x))
    oneDc <- makeOneDc(x, snapshot)
    colPerDc <- length(oneDc)
    nDc * colPerDc
}

makeOneDc <- function(x, snapshot) {
    if(x[1] %% snapshot != 0) {
        warning('Recording "on" interval must be multiple of snapshot length')
        return(NULL)
    }
    ons <- rep(snapshot, x[1]/snapshot)
    offs <- rep(snapshot, x[2] %/% snapshot)
    if(x[2] %% snapshot != 0) {
        offs <- c(offs, x[2] %% snapshot)
    }
    c(ons, offs)
}

dcToColMap <- function(x, snapshot, ncol) {
    oneDc <- makeOneDc(x, snapshot)
    rep(oneDc, length.out=ncol)
}

markOffCycle <- function(dutyCycle='2/8', snapshot=2, ncol, maxTime, shift=0) {
    # maxIx <- ceiling(maxTime/snapshot)
    dcSplit <- strsplit(dutyCycle, '/')[[1]]
    on <- as.numeric(dcSplit[1])
    off <- as.numeric(dcSplit[2])
    # if(on %% snapshot != 0 ||
    #    off %% snapshot != 0) {
    #     warning('Duty cycles are not multiples of desired snapshot length.')
    #     return(rep(FALSE, ncol)) # mark all off for this
    # }
    oneCyc <- c(
        rep(FALSE, on/snapshot),
        rep(TRUE, off %/% snapshot)
    )
    if(off %% snapshot != 0) {
        oneCyc <- c(oneCyc, TRUE)
    }
    if(shift > 0 &&
       shift < on/snapshot) {
        oneCyc <- c(oneCyc[(shift+1):length(oneCyc)], oneCyc[1:shift])
    }
    # cycleMark <- rep(
    #     c(rep(FALSE, on/snapshot), # ons marked as 1s
    #       rep(TRUE, off %/% snapshot)), # offs marked as -1s
    #     length.out=ncol
    # )
    # which(cycleMark == -1)
    numCyc <- dcToColMap(c(on, off), snapshot, ncol)
    sumCyc <- cumsum(numCyc)
    lastOn <- max(which(sumCyc < maxTime)) + 1
    markOff <- rep(oneCyc, length.out=ncol)
    if(lastOn + 1 <= ncol) {
        markOff[(lastOn+1):ncol] <- TRUE
    }
    markOff[sumCyc > maxTime] <- TRUE
    markOff
}

addRecorderInfo <- function(x, info, pascal=FALSE) {
    if(is.character(info)) {
        info <- read.csv(info, stringsAsFactors = FALSE)
    }
    if(pascal) {
        info <- tribble(
            ~station, ~dutyCycle, ~recorder,
            1, '2/8', 'ST4300',
            2, '2/8', 'ST4300',
            3, '2/8', 'ST4300',
            4, '2/8', 'ST4300',
            5, '2/8', 'ST4300',
            6, '2/0', 'ST4300',
            7, '2/0', 'SM3M',
            8, '2/2', 'SM2Bat',
            9, '2/8', 'ST4300',
            10, '2/0', 'SM3M',
            11, '2/2', 'SM2Bat',
            12, '2/8', 'ST4300',
            13, '2/0', 'SM3M',
            14, '2/2', 'SM2Bat',
            15, '2/0', 'ST4300',
            16, '2/8', 'ST4300',
            17, '2/0', 'SM3M',
            18, '2/2', 'SM2Bat',
            19, '2/8', 'ST4300',
            20, '2/8', 'ST4300',
            21, '2/8', 'ST4300',
            22, '2/8', 'ST4300'
        )
    }
    if(!all(c('station', 'dutyCycle', 'recorder') %in% colnames(info))) {
        stop('Recorder info must have columns "station", "dutyCycle" (minutes on/minutes off), and "recorder"')
    }
    hasSt <- x$station %in% info$station
    if(any(!hasSt)) {
        warning('Stations ', paste0(unique(x$station[!hasSt]), collapse=', '),
                ' are not in recorder info, names must match exactly.')
        x <- x[hasSt, ]
    }
    if(is.character(x$station) &&
       !is.character(info$station)) {
        info$station <- as.character(info$station)
    }
    x <- left_join(x, info, by='station')
    x
}

# NEED CHECK IF INTERVAL IS TIME BETWEEN STARTS OR TIME OFF CURRENTLY MIGHT BE
# not used generally - pilot project only
# formatDeployDetails <- function(x, format=c('%m/%d/%Y %H:%M:%OS', '%m-%d-%Y %H:%M:%OS',
#                                             '%Y/%m/%d %H:%M:%OS', '%Y-%m-%d %H:%M:%OS'),
#                                 sheet='deployDetails') {
#     if(is.character(x)) {
#         if(grepl('csv$', x)) {
#             x <- read.csv(x, stringsAsFactors = FALSE)
#             for(t in c('Data_Start', 'Data_End')) {
#                 x[[t]] <- gsub(' [[:alpha:]]+', '', x[[t]])
#                 x[[t]] <- lubridate::parse_date_time(x[[t]], orders = format, tz='UTC', exact=TRUE, truncated = 2)
#             }
#         } else if(grepl('xlsx$', x)) {
#             x <- read_xlsx(x, sheet=sheet, col_types='list')
#             x$Project <- as.character(x$Project)
#             x$Type <- as.character(x$Type)
#             for(t in c('Data_Start', 'Data_End')) {
#                 x[[t]] <- suppressWarnings(as.POSIXct(as.numeric(unlist(x[[t]])), origin='1970-01-01 00:00:00', tz='UTC'))
#             }
#         }
#     }
#     x <- select(x, all_of(c('Project', 'station'='DeploymentID', 'recorder'='Type',
#                             'deployTime'='Data_Start', 'retrTime'='Data_End',
#                             'dutyOn'='RecordingDuration_m', 'dutyOff'='RecordingInterval_m')))
#     x$station <- as.character(x$station)
#     x$dutyOn <- as.character(x$dutyOn)
#     x$dutyOff <- as.character(x$dutyOff)
#     x$dutyOff[x$dutyOn == 'Continuous'] <- '1'
#     x$dutyOn[x$dutyOn == 'Continuous'] <- '1'
#     x$dutyOff <- round(as.numeric(x$dutyOff) - as.numeric(x$dutyOn), 0)
#     if(any(x$dutyOff < 0)) {
#         whichNeg <- which(x$dutyOff < 0)
#         warning('Deployments ', paste0(x$station[whichNeg], collapse=', '), ' were calculated to have',
#                 ' negative "Off" cycles (RecordingInterval - RecordingDuration).')
#         
#     }
#     x$dutyCycle <- paste0(x$dutyOn, '/', x$dutyOff)
#     x$dutyOn <- NULL
#     x$dutyOff <- NULL
#     # for(t in c('deployTime', 'retrTime')) {
#     #     x[[t]] <- gsub(' [[:alpha:]]+', '', x[[t]])
#     #     x[[t]] <- lubridate::parse_date_time(x[[t]], orders = format, tz='UTC', exact=TRUE, truncated = 2)
#     # }
#     x
# }

tallySnapshots <- function(rec, stationPattern, recorderInfo, snapshot=2) {
    if(is.AcousticStudy(rec)) {
        rec <- files(rec)$recordings
        if(is.null(rec)) {
            stop('Must add recording files to AcousticStudy object.')
        }
    }
    if(!is.data.frame(rec)) {
        stop('Input must be AcousticStudy or dataframe')
    }
    if(!all(c('station', 'deployTime', 'retrTime') %in% colnames(recorderInfo))) {
        stop('recorderInfo must have columns "station", "deployTime", and "retrTime"')
    }
    if(!all(c(inherits(recorderInfo$deployTime, 'POSIXct'),
              inherits(recorderInfo$retrTime, 'POSIXct')))) {
        stop('"deployTime" and "retrTime" must be converted to POSIXct')
    }
    # rec$db <- gsub('\\.sqlite3$', '', basename(rec$db))
    rec$db <- basename(rec$db)
    # browser()
    # rec <- left_join(rec, distinct(eventSummary[c('station', 'db')]), by='db')
    rec$station <- gsub(stationPattern, '\\1', rec$db)
    naStation <- is.na(rec$station)
    if(any(naStation)) {
        nomatchDb <- unique(rec$db[naStation])
        rec <- rec[!naStation, ]
        # nomatchSta <- unique(eventSummary$station)[!(unique(eventSummary$station) %in% unique(rec$station))]
        warning('Not able to match station names to db names for databases ',
                paste0(nomatchDb, collapse=', '), '\n',
                # 'Stations ', paste0(nomatchSta, collpase=', '), ' were not paired with wav file data.')
                'These were not paired with wav file data.')
    }
    # filter by start/end
    rec <- bind_rows(lapply(split(rec, rec$station), function(s) {
        thisRec <- recorderInfo[recorderInfo$station == s$station[1], ]
        goodTime <- s$start >= thisRec$deployTime &
            s$start <= thisRec$retrTime
        if(!any(goodTime)) {
            warning('Station ', s$station[1], ' has no files between deploy/retrieval.')
        }
        s[goodTime, ]
    }))
    # rec$interval <- interval(rec$start, rec$end)
    # snapList <- vector('list', length=nrow(rec))
    # pb <- txtProgressBar(min=0, max=length(snapList))
    # for(i in seq_along(snapList)) {
    #     this <- makeTimeEvents(bin=rec[i, ], length=snapshot, units='mins',
    #                            plot=FALSE, tryFix=FALSE, progress=FALSE)$timeRange
    #     this$overlap <- this$overlap / 60
    #     this$station <- rec$station[i]
    #     # this$id <- paste0(timeEv$id[i], '_', 1:nSubs)
    #     snapList[[i]] <- this
    #     setTxtProgressBar(pb, value=i)
    # }
    #
    # snapList <- bind_rows(snapList)
    group_by(rec, station) %>%
        summarise(nSnaps = sum(length) / (60*snapshot))
    # table(rec$station)
}

makeTimeRanges <- function(start, end=20, length='2/2', units=c('secs', 'mins', 'hours')) {
    if(is.character(length)) {
        splitLen <- strsplit(length, '/')[[1]]
        onLen <- as.numeric(splitLen[1])
        if(length(splitLen) == 1) {
            offLen <- 0
        } else {
            offLen <- as.numeric(splitLen[2])
        }
    } else {
        onLen <- length
        offLen <- 0
    }
    units <- match.arg(units)
    unitScale <- switch(units,
                        'secs' = 1,
                        'mins' = 60,
                        'hours' = 3600)
    onLen <- onLen * unitScale
    offLen <- offLen * unitScale
    if(inherits(end, 'POSIXct')) {
        totLength <- as.numeric(difftime(end, start, units='secs'))
    } else {
        totLength <- end * unitScale
    }
    startSecs <- seq(from=0, to=totLength, by=onLen + offLen)
    if(startSecs[length(startSecs)] == totLength) {
        startSecs <- startSecs[-length(startSecs)]
    }
    result <- data.frame(start = start + startSecs)
    result$end <- result$start + onLen
    result$interval <- interval(result$start, result$end)
    result
}

makeTimeEvents <- function(start=NULL, end=NULL, length, units=c('secs', 'mins', 'hours'), bin, tryFix=TRUE, plot=TRUE, progress=TRUE) {
    if(is.PAMpalSettings(bin)) {
        bin <- bin@binaries$list
    }
    if(is.character(bin)) {
        binRange <- makeBinRanges(bin, progress)
    }
    if(is.data.frame(bin)) {
        binRange <- bin
    }
    if(is.null(start)) {
        start <- min(binRange$start)
    }
    if(is.null(end)) {
        end <- max(binRange$end)
    }
    # isCont <- median(as.numeric(difftime(
    #     binRange$end, binRange$start, units='secs'
    # )))
    # isCont <- isCont < 30
    # browser()
    if(isFALSE(tryFix)) {
        timeRange <- makeTimeRanges(start, end, length=length, units=units)
    } else {
        timeRange <- fixTimeRange(start=start, end=end, bin=binRange, length=length, units=units, progress=progress)
    }
    binFilt <- filter(binRange,
                      end >= timeRange$start[1],
                      start <= timeRange$end[nrow(timeRange)])
    binFilt <- checkOverlap(binFilt, timeRange)
    timeRange <- checkOverlap(timeRange, binFilt)
    if(plot) {
        nOut <- nrow(binRange) - nrow(binFilt)
        n99 <- sum(binFilt$overlapPct < .99)
        tDiff <- as.numeric(difftime(timeRange$start[2:nrow(timeRange)], timeRange$end[1:(nrow(timeRange)-1)], units=units))
        # op <- par(mfrow=c(1,3))
        # on.exit(par(op))
        # hist(binFilt$overlapPct, breaks=seq(from=0, to=1, by=.02),
        #      xlim=c(0,1),
        #      main=paste0('Pct of each binary file in event ',
        #                  '\n(', nOut, ' files outside of time range, ',
        #                  n99, ' files < .99)',
        #                  '\n(All 1 unless duty cycle mismatch btwn event/recorder)'))
        g1 <- ggplot(binFilt, aes(x=overlapPct)) +
            geom_histogram(binwidth=.02) +
            xlim(-.03, 1.03) +
            labs(title=paste0('Pct of each binary file in event ',
                              '\n(', nOut, ' files outside of time range, ',
                              n99, ' files < .99)',
                              '\n(All 1 unless duty cycle mismatch btwn event/recorder)')) +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        # hist(timeRange$overlapPct, breaks=seq(from=0, to=1, by=.02), xlim=c(0,1),
        #      main='Pct of each event with binary data\n(Should all be 1)')
        g2 <- ggplot(timeRange, aes(x=overlapPct)) +
            geom_histogram(binwidth=.02) +
            xlim(-.03, 1.03) +
            labs(title='Pct of each event with binary data\n(Should all be 1)') +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        # if(max(tDiff) < 1) {
        #     breaks <- 'Sturges'
        # } else {
        #     breaks <- 0:ceiling(max(tDiff))
        # }
        # hist(tDiff, breaks=breaks, main=paste0('Time between events (', units, ')'), xlim=c(0, max(tDiff) + 1))
        g3 <- ggplot(data.frame(timeDiff=tDiff), aes(x=timeDiff)) +
            geom_histogram(bins=50) +
            labs(title=paste0('Time between events (', units, ')')) +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        print(g1/g2/g3)
    }
    list(timeRange=timeRange, binRange=binFilt, allBin=binRange)
}

fixTimeRange <- function(start=NULL, end=NULL, bin, length='2/8', units='mins', progress=TRUE) {
    if(is.null(start)) {
        start <- min(bin$start)
    }
    if(is.null(end)) {
        end <- max(bin$end)
    }
    bin <- bin[(bin$end >= start) & (bin$start <= end), ]
    time <- makeTimeRanges(start=start, end=end, length=length, units=units)
    time <- checkOverlap(time, bin, early=TRUE)

    result <- list()
    newBin <- bin
    newTime <- time
    if(progress) {
        pb <- txtProgressBar(min=0, max=nrow(newTime), style=3)
    }
    for(i in 1:nrow(time)) {
        # cat(paste0(nrow(newTime), ' rows remaining...\n'))
        if(nrow(newTime) <= 1) {
            result[[i]] <- newTime
            if(progress) {
                setTxtProgressBar(pb, value=nrow(time))
            }
            break
        }
        # browser()
        change <- which(newTime$overlap[2:nrow(newTime)] != newTime$overlap[1:(nrow(newTime)-1)])
        if(length(change) == 0) {
            result[[i]] <- newTime
            if(progress) {
                setTxtProgressBar(pb, value=nrow(time))
            }
            break
        }
        change <- min(change) #+ 1
        result[[i]] <- newTime[1:change, ]
        lastEnd <- newTime$end[change+1]
        lastStart <- newTime$start[change+1]
        endIn <- lastEnd %within% newBin$interval &
            lastEnd != newBin$start
        startIn <- lastStart %within% newBin$interval &
            lastStart != newBin$end
        if(any(startIn)) {
            startIx <- min(which(startIn))
        } else if(any(endIn)) {
            startIx <- min(which(endIn))
        } else {
            startIx <- min(which(newBin$start >= lastEnd))
        }
        newBin <- newBin[startIx:nrow(newBin), ]
        if(nrow(newBin) == 0) {
            break
        }
        nextStart <- min(newBin$start)

        nextStart <- max(nextStart, newTime$end[change])
        # nextBinTime <- min(newBin$start)
        newTime <- makeTimeRanges(start=nextStart, end=max(newBin$end), length=length, units=units)
        newTime <- checkOverlap(newTime, newBin, early=TRUE)
        if(i/nrow(time) > .79) {
            # browser()
        }
        if(progress) {
            setTxtProgressBar(pb, value = nrow(time) - nrow(newTime))
        }
    }
    bind_rows(result)
}

checkOverlap <- function(x, y, early=FALSE) {
    if(early) {

        overlap <- rep(0, nrow(x))
        for(i in seq_along(overlap)) {
            # isOverlap <- int_overlaps(x$interval[i], y$interval)
            isOverlap <- fast_intlap(x$interval[i], y$interval)
            if(!any(isOverlap)) {
                overlap[i] <- 0
                next
            }
            # overlap[i] <- sum(int_length(intersect(x$interval[i], y$interval[isOverlap])), na.rm=TRUE)
            overlap[i] <- sum(fast_intlen(x$interval[i], y$interval[isOverlap]))
            if(i > 1 &&
               overlap[i] != overlap[i-1]) {
                break
            }
        }
        x$overlap <- overlap
    } else {
        x$overlap <- sapply(x$interval, function(i) {
            # isOverlap <- int_overlaps(i, y$interval)
            isOverlap <- fast_intlap(i, y$interval)
            if(!any(isOverlap)) {
                return(0)
            }
            # sum(int_length(intersect(i, y$interval[isOverlap])), na.rm=TRUE)
            sum(fast_intlen(i, y$interval[isOverlap]))
        })
    }
    x$overlapPct <- x$overlap / int_length(x$interval)
    x
}

fast_intlap <- function(int1, int2) {
    int1@start <= int2@start + int2@.Data & int2@start <= int1@start +
        int1@.Data
}

fast_intlen <- function(int1, int2) {
    starts <- pmax(int1@start, int2@start)
    ends <- pmin(int1@start + int1@.Data, int2@start + int2@.Data)
    spans <- as.numeric(ends) - as.numeric(starts)
    spans
}

plotOneDetFun <- function(param, model=c('C_HN', 'HN', 'HR'), add=FALSE,
                          col='black', lwd=1, title=NULL, truncDist=4e3) {
    model <- match.arg(model)
    if(model == 'HN') {
        param[2] <- 0
    }
    ranges <- seq(from=0, to=truncDist, by=10)
    if(model == 'HR') {
        prob <- 1 - exp(-(ranges/param[1])^(-param[2]))
    } else {
        prob <- 1 - (1-exp(-.5*(ranges/param[1])^2))^(exp(param[2]))
    }
    if(add) {
        lines(x=ranges, y=prob, col=col, lwd=lwd)
    } else {
        plot(x=ranges, y=prob, type='l', col=col, lwd=lwd, ylim=c(0,1.1))
        if(is.null(title)) {
            title(paste0('Param1: ', round(param[1], 0), ' Param2: ', round(param[2], 2)))
        } else {
            title(title)
        }
    }
}

plotDetFun <- function(detFun, lwd=1, truncDist=4e3, title=NULL) {
    colors <- c('black', 'red', 'blue', 'darkgreen', 'purple')
    for(i in seq_along(detFun)) {
        params <- detFun[[i]]$maxLikeParam
        
        model <- detFun[[i]]$model
        if(is.null(model)) {
            model <- ifelse(length(params) == 2, 'C_HN', 'HN')
        }
        plotOneDetFun(params, model, add= i>1, col=colors[i], lwd=lwd, title=title, truncDist=truncDist)
    }
}