#' @title
#' Estimate stemTOCvitro score
#'
#' @aliases stemTOCvitro
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return the stemTOCvitro score.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' The function will return the 0.95 upper quantile of the 629 stemTOCvitro CpGs.
#' @return The stemTOCvitro score of each sample.
#'
#' @seealso [EpiMitClocks::EpiMitClocks()]
#'
#' @references
#'
#' @examples
#' data("LungInv")
#' stemTOCvitro.v<-stemTOCvitro(data.m = bmiq.m)
#'
#' @export
#'

stemTOCvitro <- function(data.m){
  data('cugpmitclockCpG')
  common.v <- intersect(rownames(data.m),cugpmitclockCpG.v);
  print(paste("Number of represented stemTOCvitro CpGs (max=629)=",length(common.v),sep=""))
  stemTOCvitro.v <- apply(data.m[match(cugpmitclockCpG.v,rownames(data.m)),],2,quantile,0.95,na.rm=T)
  return(stemTOCvitro.v);
}

