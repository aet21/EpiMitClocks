#' @title
#' Estimate epigenetic mitotic age
#'
#' @aliases EpiMitClocks
#'
#' @description
#' This function takes as input an Illumina 450k/EPIC DNAm beta matrix and an age vector (optional) and will return the mitotic age estimates from stemTOC, stemTOCvitro, epiTOC1, epiTOC2, epiTOC3, RepliTali, HypoClock, EpiCMIT_hyper, EpiCMIT_hypo.
#'
#' @param data.m
#' DNAm beta value matrix with rows labeling Illumina 450k/EPIC CpGs and columns labeling samples.
#'
#' @param ages.v
#' Optional argument representing the chronological ages or surrogates thereof of the samples. Vector must be of same length as the number of samples.
#'
#' @details
#' The EpiMitClocks R-package implements the most widely used mitotic clocks. Some are based on age-associated hypermethylated CpGs (stemTOC, stemTOCvitro, epiTOC3, epiTOC2, epiTOC1,EpiCMIT_hyper), others based on hypomethylated CpGs (HypoClock, RepliTali, EpiCMIT_hypo).
#'
#' * epiTOC1: The epiTOC1 score is calculated as the average beta value of 385 CpGs described in Yang et al. 2014. These CpGs are promoter CpG sites that localize to Polycomb group target genes that are unmethylated in 11 different fetal tissue types.
#'
#' * epiTOC2: Building upon a dynamic model of DNA methylation gain in unmethylated CpG-rich regions, epiTOC2 can directly estimate the cumulative number of stem cell divisions in a tissue. The details of the algorithm are described in Teschendorff 2020.
#'
#' * epiTOC3: Building upon a dynamic model of DNA methylation gain in 170 unmethylated population doubling associated CpGs, epiTOC3 can directly estimate the cumulative number of stem cell divisions in a tissue.
#'
#' * stemTOCvitro: The stemTOCvitro score is calculated as the 0.95 upper quantile of the 629 stemTOCvitro CpGs. The stemTOCvitro CpGs are promoters CpGs that are unmethylated in fetal tissue-types and undergo DNA hypermethylation with increased population-doublings.
#'
#' * stemTOC: The stemTOC score is calculated as the 0.95 upper quantile of the 371 stemTOC CpGs. Compared to stemTOCvitro CpGs, the stemTOC CpGs are filtered for significant DNA hypermethylation with chronological age in large in-vivo datasets.
#'
#' * RepliTali: The RepliTali model is described in Endicott et al. 2022. It is based on 87 population doubling associated hypomethylated CpGs.
#'
#' * HypoClock: The HypoClock is based on hypomethylation at solo-WCGW sites (Yang et al. 2014). The score is calculated as 1-meanbeta(678 solo-WCGWs).
#'
#' * EpiCMIT_hyper: EpiCMIT_hyper score is calcuated as the average beta value of 184 age-associated hypermethylated CpGs (Duran-Ferrer et al. 2020).
#'
#' * EpiCMIT_hypo: EpiCMIT_hypo score is calculated as the average beta value of 1164 age-associated hypomethylated CpGs (Duran-Ferrer et al. 2020).
#'
#'
#' @return A list containing the following entries
#'
#' * epiTOC1: The epiTOC1 score of each sample.
#'
#' * epiTOC2: A list containing the following entries.
#'
#'   * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC2 model.
#'
#'   * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC2 which assumes all epiTOC2 CpGs have beta-values exactly 0 in the fetal stage.
#'
#'   * irS: This is returned only if the ages are provided, and gives the estimated average lifetime intrinsic rate of stem-cell division per sample, as derived from epiTOC2
#'
#'   * irS2: As irS, but for the approximation.
#'
#'   * irT: The median estimate over all irS values, yielding a median estimate for the intrinsic rate of stem-cell division for the tissue.
#'
#'   * irT2: As irT, but for the approximation.
#'
#' * epiTOC3: A list containing the following entries.
#'
#'   * tnsc: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using the full epiTOC3 model.
#'
#'   * tnsc2: The estimated cumulative number of stem-cell divisions per stem-cell per year and per sample using an approximation of epiTOC3 which assumes all epiTOC3 CpGs have beta-values exactly 0 in the fetal stage.
#'
#'   * irS: This is returned only if the ages are provided, and gives the estimated average lifetime intrinsic rate of stem-cell division per sample, as derived from epiTOC3.
#'
#'   * irS2: As irS, but for the approximation.
#'
#'   * irT: The median estimate over all irS values, yielding a median estimate for the intrinsic rate of stem-cell division for the tissue.
#'
#'   * irT2: As irT, but for the approximation.
#'
#'   * avETOC3: The simple average over the 170 epiTOC3 sites.
#'
#' * stemTOCvitro: The stemTOCvitro score of each sample.
#'
#' * stemTOC: The stemTOC score of each sample.
#'
#' * RepliTali: The RepliTali score of each sample.
#'
#' * HypoClock: The HypoClock score of each sample.
#'
#' * EpiCMIT_hyper: The EpiCMIT_hyper score of each sample.
#'
#' * EpiCMIT_hypo: The EpiCMIT_hypo score of each sample.
#'
#' @seealso [EpiMitClocks::epiTOC1()], [EpiMitClocks::epiTOC1()], [EpiMitClocks::epiTOC2()], [EpiMitClocks::epiTOC3()], [EpiMitClocks::RepliTali()], [EpiMitClocks::HypoClock()], [EpiMitClocks::stemTOCvitro()], [EpiMitClocks::stemTOC()], [EpiMitClocks::EpiCMIT()]
#'
#' @references
#' Yang Z, Wong A, Kuh D, et al.
#' Correlation of an epigenetic mitotic clock with cancer risk.
#' \emph{Genome Biol.} 2016
#'
#' Teschendorff AE.
#' A comparison of epigenetic mitotic-like clocks for cancer risk prediction.
#' \emph{Genome Med.} 2020
#'
#' Endicott JL, Nolte PA, Shen H, Laird PW.
#' Cell division drives DNA methylation loss in late-replicating domains in primary human cells.
#' \emph{Nat Commun.} 2022
#'
#' Duran-Ferrer M, Clot G, Nadeu F, et al.
#' The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome.
#' \emph{Nat Cancer} 2020
#'
#' @examples
#' data("LungInv")
#' EpiMitClocks.o<-EpiMitClocks(data.m = bmiq.m,ages.v = df$Age)
#'
#' @export
#'


EpiMitClocks<-function(data.m,ages.v=NULL){
  epitoc1.v<-epiTOC1(data.m)
  epitoc2.o<-epiTOC2(data.m,ages.v)
  epitoc3.o<-epiTOC3(data.m,ages.v)
  stemTOCvitro.v<-stemTOCvitro(data.m)
  stemTOC.v<-stemTOC(data.m)
  hypoclock.v<-HypoClock(data.m)
  replitali.v<-RepliTali(data.m)
  epicmit.o<-EpiCMIT(data.m)
  return(list(epiTOC1=epitoc1.v,
              epiTOC2=epitoc2.o,
              epiTOC3=epitoc3.o,
              stemTOCvitro=stemTOCvitro.v,
              stemTOC=stemTOC.v,
              RepliTali=replitali.v,
              HypoClock=hypoclock.v,
              EpiCMIT_Hyper=epicmit.o$hyperSC,
              EpiCMIT_Hypo=epicmit.o$hypoSC))
}
