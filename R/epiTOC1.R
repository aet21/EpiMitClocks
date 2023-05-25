#' @title
#' Estimate epiTOC1 score
#'
#' @aliases epiTOC1
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return the epiTOC1 score.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' The epiTOC1 score is calculated as the average beta value of 385 CpGs described in Yang et al. 2014. These CpGs are promoter CpG sites that localize to Polycomb group target genes that are unmethylated in 11 different fetal tissue types.
#'
#' @return The epiTOC1 score of each sample.
#' @seealso [EpiMitClocks::EpiMitClocks()]
#'
#' @references
#' Yang Z, Wong A, Kuh D, et al.
#' Correlation of an epigenetic mitotic clock with cancer risk.
#' \emph{Genome Biol.} 2016
#'
#'
#' @examples
#' data("LungInv")
#' epitoc1.v<-epiTOC1(data.m = bmiq.m)
#'
#' @export
#'


epiTOC1 <- function(data.m){
  data('dataETOC3')
  cpgETOC.v <- dataETOC3.l[[3]];
  common.v <- intersect(rownames(data.m),cpgETOC.v);
  print(paste("Number of represented epiTOC1 CpGs (max=385)=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
  return(pcgtAge.v);
}

