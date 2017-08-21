### merge any two TRA dataframes with colnames chr1, pos1, chr2, pos2
methodsMergeTRA <- function(dfBD="", dfDelly="", windowSize=300) 
{
  tmpBD <- dfBD
  tmpDelly <- dfDelly
  dfTRA <- rbind(tmpBD, tmpDelly)
  
  if (!is.null(dfTRA)) {
    tmpDf <- dfTRA
    names(tmpDf) <- c("chr1", "pos1", "chr2", "pos2")
    range1 <- GRanges(seqnames=tmpDf$chr1, 
                      ranges=IRanges(start=(tmpDf$pos1-windowSize/2), 
                                     end=(tmpDf$pos1+windowSize/2)))
    range2 <- GRanges(seqnames=tmpDf$chr2, 
                      ranges=IRanges(start=(tmpDf$pos2-windowSize/2), 
                                     end=(tmpDf$pos2+windowSize/2)))
    rangeReduce1 <- findOverlaps(range1, reduce(range1))
    tmpDf$clu1 <- subjectHits(rangeReduce1)
    rangeReduce2 <- findOverlaps(range2, reduce(range2))
    tmpDf$clu2 <- subjectHits(rangeReduce2)
    
    tmpDfMerge <- ddply(tmpDf, ("clu1"), function(x){
      return(ddply(x, ("clu2"), function(y){
        if (nrow(y)==1) {
          return(NULL)
        } 
        else {
          tmpDfMerging <- c(y$chr1[1], min(y$pos1), y$chr2[1], max(y$pos2), y$clu1[1], y$clu2[1])
          return(tmpDfMerging)
        }
      }))
    })
  } 
  else {
    tmpDfMerge <- NULL
  }
  
  dfRes <- na.omit(tmpDfMerge[,3:6])
  names(dfRes) <- c("chr1", "pos1", "chr2", "pos2")
  dfRes <- dfRes[order(as.numeric(dfRes$chr1), as.numeric(dfRes$pos1)),]
  attributes(dfRes) <- c(attributes(dfRes), 
                           list(method="BreakDancer, Delly")) ### may change it with respective to SV callers
  return(dfRes)
}
