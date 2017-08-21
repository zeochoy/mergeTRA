readBreakDancerTRA <- function(file="", scoreCutoff=60, readCutoff=10, windowSize=300)
{
  colClass <- c("character", "numeric", "NULL", "character", "numeric", 
                  "NULL", "character", "numeric", "numeric", "numeric")
  rawBD <- read.table(file, colClasses=colClass, fill=TRUE)
  rawBD <- rawBD[,1:8]
  names(rawBD) <- c("chr1", "pos1", "chr2", "pos2", "type", "size", 
                     "score", "reads")
  rawBD <- rawBD[rawBD$score>=scoreCutoff&
                     rawBD$reads>=readCutoff, ]
  
  ### reading in translocation
  bdTRA <- rawBD[(rawBD$type=="CTX"), ]
  
  ### merging translocation
  if (!is.null(bdTRA)) {
    tmpDf <- bdTRA[, c("chr1", "pos1", "chr2", "pos2", "reads")]
    names(tmpDf) <- c("chr1", "pos1", "chr2", "pos2", "reads")
    bdRange1 <- GRanges(seqnames=tmpDf$chr1, 
                           ranges=IRanges(start=(tmpDf$pos1-windowSize/2), 
                                          end=(tmpDf$pos1+windowSize/2)))
    bdRange2 <- GRanges(seqnames=tmpDf$chr2, 
                           ranges=IRanges(start=(tmpDf$pos2-windowSize/2), 
                                          end=(tmpDf$pos2+windowSize/2)))
    bdRangeReduce1 <- findOverlaps(bdRange1, reduce(bdRange1))
    tmpDf$clu1 <- subjectHits(bdRangeReduce1)
    bdRangeReduce2 <- findOverlaps(bdRange2, reduce(bdRange2))
    tmpDf$clu2 <- subjectHits(bdRangeReduce2)
    
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
  attributes(dfRes) <- c(attributes(dfRes), list(method="BreakDancer"))
  
  return(dfRes)
}
