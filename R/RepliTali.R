#' @title
#' Estimate RepliTali score
#'
#' @aliases RepliTali
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return the RepliTali score.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' The RepliTali model is described in Endicott et al. 2022. It is based on 87 population doubling associated hypomethylated CpGs.
#' @return The RepliTali score of each sample.
#'
#' @seealso [EpiMitClocks::EpiMitClocks()]
#'
#' @references
#' Endicott JL, Nolte PA, Shen H, Laird PW.
#' Cell division drives DNA methylation loss in late-replicating domains in primary human cells.
#' \emph{Nat Commun.} 2022
#'
#' @examples
#' data("LungInv")
#' replitali.v<-RepliTali(data.m = bmiq.m)
#' @export
#'


RepliTali<-function(data.m){
  data('Replitali')
  common.v <- intersect(rownames(data.m),replitali.cpg.v);
  print(paste("Number of represented RepliTali CpGs (max=87)=",length(common.v),sep=""))
  rep.beta.m<-data.m[match(replitali.cpg.v,rownames(data.m)),]
  rep.beta.m[1,]<-1 # Intercept term
  multi.beta.m<-rep.beta.m * replitali.coe
  replitali.v<-colSums(multi.beta.m,na.rm = T)
  return(replitali.v)
}
