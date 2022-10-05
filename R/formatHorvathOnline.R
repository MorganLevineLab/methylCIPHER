#' formatHorvathOnline
#'
#' @description A formatting Function for the Horvath Online Calculator
#'
#' @details This is a tool to prepare the Pheno and DNAm dataframes for the Horvath online calculator. This might be necessary if GrimAge estimation, or methylation-beta value based blood cell composition information is desired. Please note that while attempts have been made to ensure that the Pheno frame is reformatted properly, we caution you to check the output .csv, esp the Female and Age columns, as for any automatic formatter to ensure reasonable outputs.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#'
#' @param DNAm The methylation Beta values you will calculate epigenetic clocks with, where columns are CpGs, and rows are samples.
#' @param Pheno The phenotype data for samples, including a numeric column named "Age" and a numeric (1/0) column named "Female". Ensure that the Pheno frame has the same sample order as DNAm before running this.
#' @param writePath The path to the directory where you want the properly formatted DNAm and Pheno files (.csv) saved to be uploaded.
#'
#' @return Formatted .csv files for input to the online Horvath calculator at <http://dnamage.genetics.ucla.edu/home>. This will provide access to GrimAge, Beta value based blood cell composition, and a few epigenetic clocks.
#' @export
#'
#' @examples formatHorvathOnline(exampleBetas, examplePheno)
formatHorvathOnline <- function(DNAm, Pheno, writePath){

  #FIRST, check that the dimensions of DNAm and Pheno are correct.
  if(dim(DNAm)[1] != dim(Pheno)[1]){
    stop("Dimensions of DNAm and Pheno inputs suggest different numbers of samples! \n Please ensure you've checked that sample order is aligned in both objects.")}



  #SECOND, check that Pheno has 'Age' and 'Female' columns, or find this information.
  if(!all(c("Age","Female") %in% colnames(Pheno))){
    warning("Age and Female columns not found. Checking for columns with similar information to format.")

    if(!"Age" %in% colnames(Pheno)){
      if(sum(grepl("age", colnames(Pheno))) > 1){ #Alternative column names containing age information multiple ambiguous options
        warning("Multiple potential columns containing Age information")}

      if(sum(grepl("age", colnames(Pheno))) == 1){ #Alternative column name containing age information found
        message("Alternative column name found for 'Age'. Fixing now.")
        Pheno$Age <- Pheno[,grepl("age", colnames(Pheno))][[1]]
        Pheno <- Pheno[,-c(grepl("age", colnames(Pheno)))]}
    }

    if(!"Female" %in% colnames(Pheno)){
      if(any(grepl("female", colnames(Pheno)))){ #Checking for improper caps of Female
        message("It appears the column name 'Female' wasn't capitalized! Fixing now.")
        Pheno$Female <- Pheno[,grepl("female", colnames(Pheno))][[1]]
        Pheno <- Pheno[,-c(grepl("female", colnames(Pheno)))]}

      if(sum(grepl("Sex|sex|Male|Gender|gender", colnames(Pheno))) > 1){ #Alternative column names containing Female Sex information multiple ambiguous options
        warning("Multiple potential columns containing 'Female' information")}

      if(sum(grepl("Sex|sex|Male|Gender|gender", colnames(Pheno))) == 1){ #Alternative column name containing Female Sex information
        message("Alternative column name found for 'Female'.")
        sexType <- 1*(any(grepl("Sex|sex|Gender|gender", colnames(Pheno)))) +
          2*(any(grepl("Male", colnames(Pheno)))) # Save the type that "Female" has been converted from
        Pheno$Female <- Pheno[,grepl("Sex|sex|Male|Gender|gender", colnames(Pheno))][[1]]
        Pheno <- Pheno[,-c(grepl("Sex|sex|Male|Gender|gender", colnames(Pheno)))]}
    }

    if(!("Age" %in% colnames(Pheno))){
      warning("'Age' column couldn't be inferred Please try renaming the appropriate column manually, or add it if missing.")}
    if(!("Female" %in% colnames(Pheno))){
      warning("'Female' column couldn't be inferred Please try renaming the appropriate column manually, or add it if missing.")}
    if(!all(c("Age","Female") %in% colnames(Pheno))){
      stop("Cannot proceed without 'Age' and 'Female' columns in Pheno input!")}
  }



  #THIRD, check that the Age and Female Columns have the appropriate data in their columns
  if("sexType" %in% ls()){ #Indicates that a conversion must occur for the Female column

    if(sexType == 1){ #tends to be characters corresp to "Sex"
      Pheno$Female <- as.character(Pheno$Female) #ensure that factor variables are converted to strings

      if(any(c("F","f","M","m") %in% Pheno$Female)){ #M/F character conversion
        tempVec <- rep(NA, length(Pheno$Female))
        tempVec[Pheno$Female %in% c("F","f")] <- 1
        tempVec[Pheno$Female %in% c("M","m")] <- 0
        Pheno$Female <- as.numeric(tempVec)
        message("Column 'Female' was converted from 'M'/'F' to 0/1.")
        } else if(any(c("female","Female","FEMALE","male","Male","MALE") %in% Pheno$Female)){ # male/female character conversion
        tempVec <- rep(NA, length(Pheno$Female))
        tempVec[Pheno$Female %in% c("female", "Female", "FEMALE")] <- 1
        tempVec[Pheno$Female %in% c("male", "Male", "MALE")] <- 0
        Pheno$Female <- as.numeric(tempVec)
        message("Column 'Female' was converted from 'male'/'female' to 0/1.")
        } else stop("Not sure how to convert the original sex column to appropriate format.
                    Typically this means there was a 'gender' or 'sex' columns with 1/0 entries.")

    }

    if(sexType == 2){ #if the original sex coding was for males
      Pheno$Female <- 1 - as.numeric(Pheno$Female)
      message("Column 'Female' was recoded from males as 1's to females as 1's")}

  }
  if(class(Pheno$Age) != "numeric"){
    Pheno$Age <- as.numeric(Pheno$Age)
    if(all(is.na(Pheno$Age))){stop("Conversion of column 'Age' to numeric failed!")}
    if(any(is.na(Pheno$Age))){warning("Some Ages in column 'Age' are missing.")}
  }

  #FOURTH, format the DNAm dataframe
  DNAm <- subsetCG(DNAm, HorvathOnlineRef$Name)
  DNAm <- as.data.frame(t(DNAm)) %>% rownames_to_column("ProbeID") #rotate DNAm to match online calculator dimensions\
  message("DNA methylation data has been successfully subset to the required probes for the online calculator.")

  #FIFTH, write the csv files for the user

  write.table(DNAm,file = paste(writePath,"DNAm_Horvath_Online_Calculator_input.csv",sep="/"), row.names=F, sep="," )
  write.table(Pheno,file = paste(writePath,"Pheno_Horvath_Online_Calculator_input.csv",sep="/"),row.names=F,sep="," )

  message(paste("The .csv files 'DNAm_Horvath_Online_Calculator_input' and 'Pheno_Horvath_Online_Calculator_input' have been generated.
          They can be found in", writePath))

}



