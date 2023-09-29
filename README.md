---
title: "Introduction to R-package EpiMitClocks"
author:
- name: "Andrew E. Teschendorff"
  affiliation: 
  - CAS Key Lab of Computational Biology, PICB, SINH
date: "2023-5-25"
package: EpiMitClocks_0.1.0
output:
  BiocStyle::html_document:
    toc_float: true
---
  


# Summary

The purpose of the `EpiMitClocks` R-package is to provide functions for estimating the mitotic age of tissues from a corresponding DNA methylation (DNAm) profile. These functions are epigenetic mitotic clocks, providing separate mitotic age estimates on scales unique to each clock. So each estimate should be interpreted as a relative mitotic age, not an absolute one, except for the epiTOC2 clock which does aim to estimate the total number of stem-cell divisions in a tissue. Current version of the R-package, estimates mitotic age according to the following mitotic clocks: epiTOC, epiTOC2, HypoClock, RepliTali, epiCMIT and stemTOC. EpiMitClocks has been installed and tested on Ubuntu running R-version 4.3.1 and Windows10 running R-version 4.1.1 .

# Installation

To install:

```r
library(devtools)
devtools::install_github("aet21/EpiMitClocks")
```
Alternatively, download the EpiMitClocks-master.zip file from the Code-link above and
install locally from within R.

# References

Correlation of an epigenetic mitotic clock with cancer risk. Yang Z, Wong A, Kuh D, Paul DS, Rakyan VK, Leslie RD, Zheng SC, Widschwendter M, Beck S, Teschendorff AE. Genome Biol. 2016 Oct 3;17(1):205. doi: 10.1186/s13059-016-1064-3.

DNA methylation loss in late-replicating domains is linked to mitotic cell division. Zhou W, Dinh HQ, Ramjan Z, Weisenberger DJ, Nicolet CM, Shen H, Laird PW, Berman BP. Nat Genet. 2018 Apr;50(4):591-602. doi: 10.1038/s41588-018-0073-4.

A comparison of epigenetic mitotic-like clocks for cancer risk prediction. Teschendorff AE. Genome Med. 2020 Jun 24;12(1):56. doi: 10.1186/s13073-020-00752-3.

Cell division drives DNA methylation loss in late-replicating domains in primary human cells. Endicott JL, Nolte PA, Shen H, Laird PW. Nat Commun. 2022 Nov 8;13(1):6659. doi: 10.1038/s41467-022-34268-8.

The proliferative history shapes the DNA methylome of B-cell tumors and predicts clinical outcome. Duran-Ferrer M, Clot G, Nadeu F, Beekman R, Baumann T, Nordlund J, Marincevic-Zuniga Y, Lönnerholm G, Rivas-Delgado A, Martín S, Ordoñez R, Castellano G, Kulis M, Queirós AC, Lee ST, Wiemels J, Royo R, Puiggrós M, Lu J, Giné E, Beà S, Jares P, Agirre X, Prosper F, López-Otín C, Puente XS, Oakes CC, Zenz T, Delgado J, López-Guillermo A, Campo E, Martín-Subero JI. Nat Cancer. 2020 Nov;1(11):1066-1081. doi: 10.1038/s43018-020-00131-2

 
