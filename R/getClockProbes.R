#' getClockProbes
#'
#' @param DNAm The methylation Beta values you will calculate epigenetic clocks with, where columns are CpGs, and rows are samples.
#'
#' @return A table which compares the number of probes required by the full clocks in the current package, and the ones in your data. If you are missing many probes for a clock, imputation is necessary, though you might consider using a different clock or dataset if possible.
#' @export
#'
#' @examples getClockProbes(exampleBetas)
getClockProbes <- function(DNAm){

  ClockDataList <- as.data.frame(data(package = "methylCIPHER")[3][[1]])$Item
  ClockDataList <- stringr::str_subset(ClockDataList, pattern = "CpG")

  totalProbes <- vector(mode = "numeric", length = length(ClockDataList))
  presentProbes <- vector(mode = "numeric", length = length(ClockDataList))
  percentPresent <- vector(mode = "character", length = length(ClockDataList))

  for(i in 1:length(ClockDataList)){

    x <- get(ClockDataList[i])

    if(ClockDataList[i] %in% c("Bocklandt_CpG","EpiToc_CpGs","hypoClock_CpGs","Garagnani_CpG","Weidner_CpGs")){

      currentCpGList <- x

    } else if(ClockDataList[i] == "EpiToc2_CpGs"){

      currentCpGList <- rownames(x)

    } else if(ClockDataList[i] %in% c("DNAmTL_CpGs", "HRSInCHPhenoAge_CpGs", "Horvath1_CpGs",
                                      "MiAge_CpGs", "PEDBE_CpGs", "Zhang2019_CpGs")){

      pull <- function(x1,y1) {x1[,if(is.name(substitute(y1))) deparse(substitute(y1)) else y1, drop = FALSE][[1]]}

      CpGColumn <- grepl("CpG|Marker|ID|id", colnames(x))
      currentCpGList <- pull(x,colnames(x)[CpGColumn])

    } else {

      CpGColumn <- grepl("CpG|Marker|ID|id", colnames(x))

      currentCpGList <- x[CpGColumn][,1]

    }

    totalProbes[i] <- length(currentCpGList)
    presentProbes[i] <- sum(currentCpGList %in% colnames(DNAm))
    percentPresent[i] <- paste(round((presentProbes[i]/totalProbes[i])*100,0),"%",sep = "")

  }

  ProbeTable <- data.frame(Clock = sapply(strsplit(ClockDataList,"_"), `[`, 1),
                           `Total Probes` = totalProbes,
                           `Present Probes` = presentProbes,
                           `Percent Present` = percentPresent)

  ProbeTable

}
