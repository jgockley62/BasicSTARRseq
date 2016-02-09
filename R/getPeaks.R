################################################################################
############## Intern function to calculate potential peak positions ###########

.getPeakPos <- function(sampleCov, quan, peakWidth){
    which.quan <- which(runValue(sampleCov) > quan)
    
    runVal <- runValue(sampleCov)[which.quan]
    potSummitPos <- (start(sampleCov)[which.quan]+
                         end(sampleCov)[which.quan])%/%2
    runValSortIx <- order(runVal, decreasing=TRUE)
    peakIx <- rep(FALSE, length(runValSortIx))
    nopeakIx <- rep(FALSE, length(runValSortIx))
    
    for (i in runValSortIx){
        if (!nopeakIx[i]){
            peakIx[i] <- TRUE
            j <- i + 1
            while (j <= length(potSummitPos) &&
                   potSummitPos[j] - potSummitPos[i] < peakWidth){
                nopeakIx[j] <- TRUE
                j <- j + 1
            }
            j <- i - 1
            while (j >= 1 &&
                   potSummitPos[i]-potSummitPos[j] < peakWidth){
                nopeakIx[j] <- TRUE
                j <- j - 1
            }
        }
    }
    return(potSummitPos[peakIx])
}

################################################################################
### Intern function to correct enrichment to limits of confidence intervals ####

.correctConfInt <- function(x, total1, total2){
    if (x[1] == x[2])
        return(1)
    pt1 <- prop.test(x=x[1]*total1, n=total1)$conf.int
    pt2 <- prop.test(x=x[2]*total2, n=total2)$conf.int
    if (x[1] > x[2]){
        return(max(pt1[1]/pt2[2], 1))
    } else {
        return(min(pt1[2]/pt2[1], 1))
    }
}

################################################################################
########## Basic peak calling on STARR-seq data object #########################

setMethod(f="getPeaks",
          signature="STARRseqData",
          definition=function(object, minQuantile=0.9, peakWidth=500,
                              maxPval=0.001, deduplicate=TRUE, model=1){
              if (minQuantile <= 0 || minQuantile > 1){
                  message("function getPeaks parameter minQuantile ",
                          "(must be between 0 and 1) set to 0.9")
                  minQuantile <- 0.9
              }
              if (peakWidth <= 0){
                  message("function getPeaks parameter peakWidth ",
                          "(must be greater than 0) set to 500")
                  peakWidth <- 500
              }
              if (maxPval < 0 || maxPval > 1){
                  message("function getPeaks parameter maxPval ",
                          "(must be between 0 and 1) set to 0.001")
                  maxPval <- 0.001
              }

              if (deduplicate) {
                message("deduplicate data")
                sampleRu <- unique(sample(object))
                controlRu <- unique(control(object))
              } else {
                  sampleRu <- sample(object)
                  controlRu <- control(object)
              }

              message("calculate coverage of data")
              sampleCov <- coverage(sampleRu)
              controlCov <- coverage(controlRu)

              cntSample <- length(sampleRu)
              cntControl <- length(controlRu)
              fraction <- cntSample/(cntControl+cntSample)

              # potential summit positions
              message("calculate potential summit positions")
              potSummitPos <- lapply(sampleCov, .getPeakPos,
                                    quantile(unlist(sampleCov), minQuantile), 
                                    peakWidth)

              summitsStarts <- pmax(1, unlist(potSummitPos)-
                                        ceiling(peakWidth/2))
              summitsStarts <- pmin(summitsStarts, 
                                    as.vector(Rle(seqlengths(sampleRu),
                      sapply(potSummitPos,length))-peakWidth+1))
              summits <- GRanges(seqnames=Rle(names(potSummitPos),
                                              sapply(potSummitPos, length)),
                               ranges=IRanges(start=summitsStarts,
                                                width=peakWidth),
                               seqinfo=seqinfo(sampleRu))

              # summits control coverage / summits sample coverage
              message(length(summits), " potential summits found,",
                            " calculate coverage of summits")
              summitSampleCov <- vector("list", length(sampleCov))
              summitControlCov <- vector("list", length(sampleCov))
              for (chr in seq_along(sampleCov)){
                  summitSampleCov[[chr]] <- as.vector(
                      sampleCov[[chr]][potSummitPos[[chr]]])
                  summitControlCov[[chr]] <- as.vector(
                      controlCov[[chr]][potSummitPos[[chr]]])
              }
              summits$sampleCov <- unlist(summitSampleCov)
              summits$controlCov <- unlist(summitControlCov)

              # approximation summits p-Val (with control coverage)
              message("calculate p-Value of summits")
              if (model == 1){
                  binModel <- data.frame(
                      nrSuccess=summits$sampleCov,
                      nrTrials=cntSample,
                      estProp=summits$controlCov/cntControl)
              } else {
                  binModel <- data.frame(
                      nrSucess=summits$sampleCov,
                      nrTrials=summits$sampleCov+summits$controlCov,
                      estProp=cntSample/(cntSample+cntControl))

              }
              summits$pVal <- apply(binModel, 1, function(y)
                  pbinom(q=y[1]-1, size=y[2], prob=y[3], lower.tail=FALSE))

              # filter p-Val
              summits <- summits[summits$pVal <= maxPval]

              # summits medianControlCov
              summitMedianControlCov <- vector("list", length(sampleCov))
              for (chr in seq_along(sampleCov)){
                  sidx <- summits[seqnames(summits) == names(sampleCov)[chr]]
                  rng <- ranges(sidx)
                  m <- matrix(controlCov[[chr]][rng], peakWidth)
                  summitMedianControlCov[[chr]] <- apply(m, 2, median)
              }
              summits$medianControlCov <- round(unlist(summitMedianControlCov))

              # corrected summits p-Val (with max(controlCov, medianControlCov))
              where.correct <- summits$medianControlCov > summits$controlCov
              if (sum(where.correct) > 0){
                  if (model == 1){
                      binModel <- data.frame(
                          nrSuccess=summits$sampleCov[where.correct],
                          nrTrials=cntSample,
                          estProp=summits$medianControlCov[where.correct]/
                              cntControl)
                  } else {
                      binModel <- data.frame(
                          nrSucess=summits$sampleCov[where.correct],
                          nrTrials=(summits$sampleCov+
                                    summits$medianControlCov)[where.correct],
                          estProp=cntSample/(cntSample+cntControl))

                  }
                  summits$pVal[where.correct] = apply(binModel, 1, function(y)
                      pbinom(q=y[1]-1, size=y[2], prob=y[3], lower.tail=FALSE))
              }

              # filter p-Val
              summits <- summits[summits$pVal <= maxPval]

              summits$controlCov <- pmax(summits$controlCov,
                                        summits$medianControlCov)

              # summits enrichment
              message(length(summits), " summits with p-Value <= ",
                            maxPval, ", calculate enrichment")
              summits$enrichment <- apply(cbind(summits$sampleCov/cntSample,
                                               summits$controlCov/cntControl),
                                         1, .correctConfInt, total1=cntSample,
                                         total2=cntControl)
              return(summits[, -4])
          }
)
