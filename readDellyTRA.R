readDellyTrx <- function(dataDir=".", readCutoff=10, scoreCutoff=30, pass=TRUE, windowSize=300){
  ### reading in SV predictions
  dellyFileList <- list.files(dataDir, full.names=TRUE)
  
  rawDelly <- lapply(dellyFileList, function(x) {
    dellyData <- try(read.table(x, as.is=TRUE), silent=TRUE)
    if (is.data.frame(dellyData)) {
      dellyData <- dellyData[, c(1:2, 7:8, 10)]
      names(dellyData) <- c("chr1", "start", "pass", 
                            "info", "detail")
      return(dellyData)
    } else {
      return(NULL)
    }
  })
  
  dfDelly <- do.call(rbind, rawDelly)
  
  if (pass) {
    dfDelly <- dfDelly[dfDelly$pass=="PASS", ]
  }
  
  dfDelly$chr2 <- as.numeric(gsub(".+;CHR2=(\\d+).*", "\\1", dfDelly$info))
  dfDelly$end <- as.numeric(gsub(".+;END=(\\d+);.+", "\\1", dfDelly$info))
  dfDelly$type <- gsub(".+;SVTYPE=([A-Z]+);.+", "\\1", dfDelly$info)
  dfDelly$reads <- as.numeric(gsub(".+;PE=(\\d+);.+", "\\1", dfDelly$info))
  dfDelly$score <- as.numeric(gsub(".+;MAPQ=(\\d+).*", "\\1", dfDelly$info))
  
  dfDelly <- dfDelly[dfDelly$score>=scoreCutoff & dfDelly$reads>=readCutoff, ]
  
  dfDelly <- dfDelly[, c("chr1", "start", "chr2", "end", "type","reads", "score")]
  
  ### reading in translocation
  dfDellyTRA <- dfDelly[dfDelly$type=="BND", ]
  
  ### merging translocation
  if (!is.null(dfDellyTRA)) {
    tmpDf <- dfDellyTRA[, c("chr1", "start", "chr2", "end", "reads")]
    names(tmpDf) <- c("chr2", "pos2", "chr1", "pos1", "reads")
    tmpDf <- tmpDf[,c("chr1", "pos1", "chr2", "pos2", "reads")]
    dellyRange1 <- GRanges(seqnames=tmpDf$chr1, 
                           ranges=IRanges(start=(tmpDf$pos1-windowSize/2), 
                                          end=(tmpDf$pos1+windowSize/2)))
    dellyRange2 <- GRanges(seqnames=tmpDf$chr2, 
                           ranges=IRanges(start=(tmpDf$pos2-windowSize/2), 
                                          end=(tmpDf$pos2+windowSize/2)))
    dellyRangeReduce1 <- findOverlaps(dellyRange1, reduce(dellyRange1))
    tmpDf$clu1 <- subjectHits(dellyRangeReduce1)
    dellyRangeReduce2 <- findOverlaps(dellyRange2, reduce(dellyRange2))
    tmpDf$clu2 <- subjectHits(dellyRangeReduce2)
    tmpDfMerge <- ddply(tmpDf, ("clu1"), function(x){
      return(ddply(x, ("clu2"), function(y){
        if (nrow(y)==1) {
          return(y)
        } 
        else {
          tmpDfMerging <- c(y$chr1[1], min(y$pos1), y$chr2[1], max(y$pos2), max(y$reads), y$clu1[1], y$clu2[1])
          return(tmpDfMerging)
          #return(y[which.max(y$rp_support),])
        }
      }))
    })
  }
  else {
    tmpDfMerge <- NULL
  }
  
  ### return results
  dfRes <- na.omit(tmpDfMerge[,c("chr1", "pos1", "chr2", "pos2")])
  dfRes <- dfRes[order(dfRes$chr1, dfRes$pos1),]
  attributes(dfRes) <- c(attributes(dfRes), list(method="DELLY"))
  
  return(dfRes)
}
