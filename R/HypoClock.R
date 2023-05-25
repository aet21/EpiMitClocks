#' @title
#' Estimate HypoClock score
#'
#' @aliases HypoClock
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return the HypoClock score.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' The HypoClock is based on hypomethylation at solo-WCGW sites. The score is calculated as 1-meanbeta(678 solo-WCGWs).
#'
#' @return The HypoClock score of each sample.
#'
#' @seealso [EpiMitClocks::EpiMitClocks()]
#'
#' @references
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' @examples
#' data("LungInv")
#' hypoclock.v<-HypoClock(data.m = bmiq.m)
#'
#' @export
#'

HypoClock <- function(data.m){
  data('dataETOC3')
  soloCpG.v <- dataETOC3.l[[4]];
  common.v <- intersect(rownames(data.m),soloCpG.v);
  print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  hypoSC.v <- 1-colMeans(data.m[map.idx,],na.rm=TRUE);
  return(hypoSC.v);
}
