#' @title
#' Estimate epiTOC2 scores
#'
#' @aliases epiTOC2
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an age vector (optional) and will return the epiTOC2 scores.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @param ages.v
#' Optional argument representing the chronological ages or surrogates thereof of the samples. Vector must be of same length as the number of samples.
#' @details
#' Building upon a dynamic model of DNA methylation gain in unmethylated CpG-rich regions, epiTOC2 can directly estimate the cumulative number of stem cell divisions in a tissue. The details of the algorithm are described in Teschendorff et al. 2020.
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC2 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC2 which assumes all epiTOC2 CpGs have beta-values exactly 0 in the fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the estimated average lifetime intrinsic rate of stem-cell division per sample, as derived from epiTOC2
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
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
#' epitoc2.o<-epiTOC2(data.m = bmiq.m,ages.v = df$Age)
#'
#' @export
#'

epiTOC2<-function(data.m,ages.v=NULL){
  data('dataETOC3')
  estETOC2.m <- dataETOC3.l[[2]]
  map.idx <- match(rownames(estETOC2.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
  tmp.m <- data.m[map.idx[rep.idx],];
  TNSC.v <- 2*colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
  TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
  estIR.v <- NULL; estIR2.v <- NULL;
  estIR <- NULL;  estIR2 <- NULL;
  if(!is.null(ages.v)){
    estIR.v <- TNSC.v/ages.v;
    estIR <- median(estIR.v,na.rm=TRUE);
    estIR2.v <- TNSC2.v/ages.v;
    estIR2 <- median(estIR2.v,na.rm=TRUE);
  }
  return(list(tnsc=TNSC.v,tnsc2=TNSC2.v,irS=estIR.v,irS2=estIR2.v,irT=estIR,irT2=estIR2));

}
