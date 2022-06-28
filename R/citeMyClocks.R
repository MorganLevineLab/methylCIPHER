#' citeMyClocks
#'
#' @param clockList A character value or vector naming the function(s) you used to calculate epigenetic clocks.
#' @param prettyprint If you want the console output to look nice. Default is true. If you would like to store citations as a vector for export or future reference, then make FALSE.
#'
#' @return This will return a vector list of the appropriate citations printed to your console, or you can save the vector for later reference.
#' @export
#'
#' @examples citeMyClocks("calcPhenoAge")
citeMyClocks <- function(clockList, prettyprint = TRUE){

  if(class(clockList) != "character"){

    stop("Please ensure that your clockList variable is a character string or character vector")

  }

  citations <- "Clock Citations:"

  if("calcCoreClocks" %in% clockList){

    citations <- c(citations, "Hannum, G., Guinney, J., Zhao, L., Zhang, L., Hughes, G., Sadda, S. V., … Zhang, K. (2013). \n Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. \n Molecular Cell, 49(2), 359–367. https://doi.org/10.1016/j.molcel.2012.10.016",
                   "Horvath, S. (2013). \n DNA methylation age of human tissues and cell types. \n Genome Biology, 14(10), 3156. https://doi.org/10.1186/gb-2013-14-10-r115",
                   "Levine, M. E., Lu, A. T., Quach, A., Chen, B. H., Assimes, T. L., Bandinelli, S., … Horvath, S. (2018). \n An epigenetic biomarker of aging for lifespan and healthspan. \n Aging, 10(4), 573–591. https://doi.org/10.18632/aging.101414",
                   "Teschendorff, A. E. (2020). \n A comparison of epigenetic mitotic-like clocks for cancer risk prediction. \n Genome Medicine, 12(1), 56. https://doi.org/10.1186/s13073-020-00752-3")}

  if(any(c("calcDunedinPoAm38", "DunedinPoAm38::PoAmProjector","DunedinPoAm38::PoAmProjector()") %in% clockList)){
    citations <- c(citations, "Belsky, D. W., Caspi, A., Arseneault, L., Baccarelli, A., Corcoran, D. L., Gao, X., … Moffitt, T. E. (2020). \n Quantification of the pace of biological aging in humans through a blood test, the DunedinPoAm DNA methylation algorithm. \n ELife, 9, e54870. https://doi.org/10.7554/eLife.54870")}

  if(any(c("calcDunedinPACE", "DunedinPACE::PoAmProjector","DunedinPACE::PoAmProjector()") %in% clockList)){
    citations <- c(citations, "Belsky, D. W., Caspi, A., Corcoran, D. L., Sugden, K., Poulton, R., Arseneault, L., … Moffitt, T. E. (2022). \n DunedinPACE, a DNA methylation biomarker of the pace of aging.\n  ELife, 11, e73420. https://doi.org/10.7554/eLife.73420")}

  if("calcBocklandt" %in% clockList){
    citations <- c(citations, "Bocklandt, S., Lin, W., Sehl, M. E., Sánchez, F. J., Sinsheimer, J. S., Horvath, S., & Vilain, E. (2011). \n Epigenetic Predictor of Age. PLOS ONE, 6(6), e14821. \n https://doi.org/10.1371/journal.pone.0014821")}

  if("calcBohlin" %in% clockList){
    citations <- c(citations, "Bohlin, J., Håberg, S. E., Magnus, P., Reese, S. E., Gjessing, H. K., Magnus, M. C., … Nystad, W. (2016). \n Prediction of gestational age based on genome-wide differentially methylated regions. \n Genome Biology, 17(1), 1–9. https://doi.org/10.1186/s13059-016-1063-4")
  }

  if("calcGaragnani" %in% clockList){
    citations <- c(citations, "Garagnani, P., Bacalini, M. G., Pirazzini, C., Gori, D., Giuliani, C., Mari, D., … Franceschi, C. (2012). \n Methylation of ELOVL2 gene as a new epigenetic marker of age. \n Aging Cell, 11(6), 1132–1134. https://doi.org/10.1111/acel.12005")}

  if("calcHannum" %in% clockList){
    citations <- c(citations, "Hannum, G., Guinney, J., Zhao, L., Zhang, L., Hughes, G., Sadda, S. V., … Zhang, K. (2013). \n Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. \n Molecular Cell, 49(2), 359–367. https://doi.org/10.1016/j.molcel.2012.10.016")}

  if("calcHorvath1" %in% clockList){
    citations <- c(citations, "Horvath, S. (2013). \n DNA methylation age of human tissues and cell types. \n Genome Biology, 14(10), 3156. https://doi.org/10.1186/gb-2013-14-10-r115")}

  if("calcHorvath2" %in% clockList){
    citations <- c(citations, "Horvath, S., Oshima, J., Martin, G. M., Lu, A. T., Quach, A., Cohen, H., … Raj, K. (2018). \n Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies. \n Aging, 10(7), 1758–1775. https://doi.org/10.18632/aging.101508")}

  if("calcKnight" %in% clockList){
    citations <- c(citations, "Knight, A. K., Craig, J. M., Theda, C., Bækvad-Hansen, M., Bybjerg-Grauholm, J., Hansen, C. S., … Smith, A. K. (2016).\n An epigenetic clock for gestational age at birth based on blood methylation data. \n Genome Biology, 17(1), 1–11. https://doi.org/10.1186/s13059-016-1068-z")}

  if(any(c("calcLeeRobust","calcLeeControl","calcLeeRefinedRobust") %in% clockList)){
    citations <- c(citations, "Lee, Y., Choufani, S., Weksberg, R., Wilson, S. L., Yuan, V., Burt, A., … Horvath, S. (2019). \n Placental epigenetic clocks: Estimating gestational age using placental DNA methylation levels. \n Aging, 11(12), 4238–4253. https://doi.org/10.18632/aging.102049")}

  if("calcPhenoAge" %in% clockList){
    citations <- c(citations, "Levine, M. E., Lu, A. T., Quach, A., Chen, B. H., Assimes, T. L., Bandinelli, S., … Horvath, S. (2018). \n An epigenetic biomarker of aging for lifespan and healthspan. \n Aging, 10(4), 573–591. https://doi.org/10.18632/aging.101414")}

  if("calcLin" %in% clockList){
    citations <- c(citations, "Lin, Q., Weidner, C. I., Costa, I. G., Marioni, R. E., Ferreira, M. R. P., Deary, I. J., & Wagner, W. (2016). \n DNA methylation levels at individual age-associated CpG sites can be indicative for life expectancy. \n Aging, 8(2), 394–401. https://doi.org/10.18632/aging.100908")}

  if("calcDNAmTL" %in% clockList){
    citations<- c(citations, "Lu, A. T., Seeboth, A., Tsai, P. C., Sun, D., Quach, A., Reiner, A. P., Kooperberg, C., … Horvath, S. (2019). \n DNA methylation-based estimator of telomere length. \n Aging, 11(16), 5895–5923. https://doi.org/10.18632/aging.102173")}

  if("calcMayne" %in% clockList){
    citations <- c(citations, "Mayne, B. T., Leemaqz, S. Y., Smith, A. K., Breen, J., Roberts, C. T., & Bianco-Miotto, T. (2017). \n Accelerated placental aging in early onset preeclampsia pregnancies identified by DNA methylation. \n Epigenomics, 9(3), 279–289. https://doi.org/10.2217/epi-2016-0103")}

  if(any(c("calcAlcoholMcCartney","calcBMIMcCartney","calcSmokingMcCartney") %in% clockList)){
    citations <- c(citations, "McCartney, D. L., Hillary, R. F., Stevenson, A. J., Ritchie, S. J., Walker, R. M., Zhang, Q., … Marioni, R. E. (2018). \n Epigenetic prediction of complex traits and death. \n Genome Biology, 19(1), 136. https://doi.org/10.1186/s13059-018-1514-1")}

  if("calcPEDBE" %in% clockList){
    citations<- c(citations, "McEwen, L. M., O’Donnell, K. J., McGill, M. G., Edgar, R. D., Jones, M. J., MacIsaac, J. L., … Kobor, M. S. (2020). \n The PedBE clock accurately estimates DNA methylation age in pediatric buccal cells. \n PNAS, 117(38), 23329–23335. https://doi.org/10.1073/pnas.1820843116")}

  if("calcDNAmClockCortical" %in% clockList){
    citations<- c(citations,"Shireby, G. L., Davies, J. P., Francis, P. T., Burrage, J., Walker, E. M., Neilson, G. W. A., … Mill, J. (2020). \n Recalibrating the epigenetic clock: implications for assessing biological age in the human cortex. \n Brain, 1–13. https://doi.org/10.1093/brain/awaa334")}

  if(any(c("calcEpiTOC2","calcHypoClock") %in% clockList)){
    citations <- c(citations, "Teschendorff, A. E. (2020). \n A comparison of epigenetic mitotic-like clocks for cancer risk prediction. \n Genome Medicine, 12(1), 56. https://doi.org/10.1186/s13073-020-00752-3")}

  if("calcVidalBralo" %in% clockList){
    citations <- c(citations, "Vidal-Bralo, L., Lopez-Golan, Y., & Gonzalez, A. (2016). \n Simplified assay for epigenetic age estimation in whole blood of adults. \n Frontiers in Genetics, 7(JUL), 1–7. https://doi.org/10.3389/fgene.2016.00126")}

  if("calcWeidner" %in% clockList){
    citations <- c(citations, "Weidner, C. I., Lin, Q., Koch, C. M., Eisele, L., Beier, F., Ziegler, P., … Wagner, W. (2014). \n Aging of blood can be tracked by DNA methylation changes at just three CpG sites. \n Genome Biology, 15(2). https://doi.org/10.1186/gb-2014-15-2-r24")}

  if("calcEpiTOC" %in% clockList){
    citations <- c(citations, "Yang, Z., Wong, A., Kuh, D., Paul, D. S., Rakyan, V. K., Leslie, R. D., … Teschendorff, A. E. (2016). \n Correlation of an epigenetic mitotic clock with cancer risk. \n Genome Biology, 17(1), 1–18. https://doi.org/10.1186/s13059-016-1064-3")}

  if("calcMiAge" %in% clockList){
    citations <- c(citations, "Youn, A., & Wang, S. (2018). \n The MiAge Calculator: a DNA methylation-based mitotic age calculator of human tissue types. \n Epigenetics, 13(2), 192–206. https://doi.org/10.1080/15592294.2017.1389361")}

  if("calcZhang" %in% clockList){
    citations <- c(citations, "Zhang, Y., Wilson, R., Heiss, J., Breitling, L. P., Saum, K. U., Schöttker, B., … Brenner, H. (2017). \n DNA methylation signatures in peripheral blood strongly predict all-cause mortality. \n Nature Communications, 8. https://doi.org/10.1038/ncomms14617")}

  if("calcZhang2019" %in% clockList){
    citations <- c(citations, "Zhang, Q., Vallerga, C. L., Walker, R. M., Lin, T., Henders, A. K., Montgomery, G. W., … Visscher, P. M. (2019). \n Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing. \n Genome Medicine, 11(1), 54. https://doi.org/10.1186/s13073-019-0667-1")}

  if("HRSInChPhenoAge" %in% clockList){
    citations <- c(citations, "HRSInChPhenoAge Citation Forthcoming")}

  if(any(c("prcPhenoAge::calcPRCPhenoAge","prcPhenoAge::calcnonPRCPhenoAge", "calcPRCPhenoAge", "calcnonPRCPhenoAge") %in% clockList)){
    citations <- c(citations, "prcPhenoAge clocks Citation Forthcoming")}

  if(prettyprint){

    writeLines(citations)

  } else print(citations)


}
