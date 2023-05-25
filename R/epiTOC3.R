#' @title
#' Estimate epiTOC3 scores
#'
#' @aliases epiTOC3
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an age vector (optional) and will return the epiTOC3 scores.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @param ages.v
#' Optional argument representing the chronological ages or surrogates thereof of the samples. Vector must be of same length as the number of samples.
#'
#' @details
#' Building upon a dynamic model of DNA methylation gain in 170 unmethylated population doubling associated CpGs, epiTOC3 can directly estimate the cumulative number of stem cell divisions in a tissue.
#'
#' @return A list containing the following entries
#'
#' * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC3 model.
#' * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC3 which assumes all epiTOC3 CpGs have beta-values exactly 0 in the fetal stage.
#' * irS: This is returned only if the ages are provided, and gives the estimated average lifetime intrinsic rate of stem-cell division per sample, as derived from epiTOC3
#' * irS2: As irS, but for the approximation.
#' * irT: The median estimate over all irS values, yielding a median estimate for the intrinsic rate of stem-cell division for the tissue.
#' * irT2: As irT, but for the approximation.
#' * avETOC3: The simple average over the 170 epiTOC3 sites.
#'
#' @seealso [EpiMitClocks::EpiMitClocks()]
#'
#' @references
#'
#'
#' @examples
#' data("LungInv")
#' epitoc3.o<-epiTOC3(data.m = bmiq.m,ages.v = df$Age)
#'
#' @export
#'


epiTOC3<-function(data.m,ages.v=NULL){
  data('dataETOC3')
  estETOC3.m <- dataETOC3.l[[1]];
  map.idx <- match(rownames(estETOC3.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Number of represented epiTOC3 CpGs (max=170)=",length(rep.idx),sep=""))
  tmp.m <- data.m[map.idx[rep.idx],];
  avETOC3.v <- colMeans(tmp.m);
  TNSC3.v <- 2*colMeans(diag(1/(estETOC3.m[rep.idx,1]*(1-estETOC3.m[rep.idx,2]))) %*% (tmp.m - estETOC3.m[rep.idx,2]),na.rm=TRUE);
  TNSC32.v <- 2*colMeans(diag(1/estETOC3.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
  estIR3.v <- NULL; estIR32.v <- NULL;
  estIR3 <- NULL;  estIR32 <- NULL;
  if(!is.null(ages.v)){
    estIR3.v <- TNSC3.v/ages.v;
    estIR3 <- median(estIR3.v,na.rm=TRUE);
    estIR32.v <- TNSC32.v/ages.v;
    estIR32 <- median(estIR32.v,na.rm=TRUE);
  }
  return(list(tnsc=TNSC3.v,tnsc2=TNSC32.v,irS=estIR3.v,irS2=estIR32.v,irT=estIR3,irT2=estIR32,avETOC3=avETOC3.v))
}
