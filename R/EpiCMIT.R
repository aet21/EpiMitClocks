#' @title
#' Estimate EpiCMIT scores
#'
#' @aliases EpiCMIT
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and will return EpiCMIT scores.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @details
#' EpiCMIT is based on two groups of CpGs, a 184 age-associated hypermethylated CpG list and a 1164 hypomethylated CpG list (Duran-Ferrer et al. 2020) . The EpiCMIT function calculates the average DNAm of the two groups separately and returns the EpiCMIT_hyper and EpiCMIT_hypo scores.
#' @return A list containing the following entries.
#'
#' * hyperSC: EpiCMIT-hyper score of each sample.
#'
#' * hypoSC: EpiCMIT-hypo score of each sample.
#'
#' @seealso [EpiMitClocks::EpiMitClocks()]
#' @references
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#'
#' @examples
#' data("LungInv")
#' epicmit.o<-EpiCMIT(data.m = bmiq.m)
#' @export
#'

EpiCMIT <- function(data.m){
  data('EpiCMITcpgs')
  epiCMIT.lv <- list();
  epiCMIT.lv[[1]] <- as.vector(epiCMIT.df[grep("hyper",epiCMIT.df[,7]),1]);
  epiCMIT.lv[[2]] <- as.vector(epiCMIT.df[grep("hypo",epiCMIT.df[,7]),1]);
  nhyper <- length(intersect(epiCMIT.lv[[1]],rownames(data.m)));
  nhypo <- length(intersect(epiCMIT.lv[[2]],rownames(data.m)));
  print(paste("Number of represented EpiCMIT_hyper CpGs (max=184)=",nhyper,sep=""))
  print(paste("Number of represented EpiCMIT_hypo CpGs (max=1164)=",nhypo,sep=""))
  hyperSC.v <- colMeans(data.m[match(intersect(epiCMIT.lv[[1]],rownames(data.m)),rownames(data.m)),]);
  hypoSC.v <- 1- colMeans(data.m[match(intersect(epiCMIT.lv[[2]],rownames(data.m)),rownames(data.m)),]);

  return(list(hyperSC=hyperSC.v,hypoSC=hypoSC.v));
}
