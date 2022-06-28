## code to prepare clock datasets goes here
library(dplyr)
library(glmnet)

#PhenoAge Prep
PhenoAge_CpGs <-as.data.frame(readr::read_csv("data-raw/Clock_Levine_513_CpG.csv"))
PhenoAge_CpGs <- PhenoAge_CpGs%>%
  select(CpG, Weight)
usethis::use_data(PhenoAge_CpGs, overwrite = TRUE)

#Horvath prep
Horvath1_CpGs <- readr::read_csv("data-raw/13059_2013_3156_MOESM23_ESM.csv")
Horvath1_CpGs <- Horvath1_CpGs %>%
  select(CpGmarker, CoefficientTraining) %>%
  filter(CpGmarker != "(Intercept)")
usethis::use_data(Horvath1_CpGs, overwrite = TRUE)

Horvath2_CpGs <- as.data.frame(readr::read_csv("data-raw/datSkinClockRecalibratedApril2018.csv"))
Horvath2_CpGs <- Horvath2_CpGs %>%
  select(ID, Coef) %>%
  filter(ID != "(Intercept)")
usethis::use_data(Horvath2_CpGs, overwrite = TRUE)

#Hannum prep
Hannum_CpGs <- as.data.frame(readxl::read_excel("data-raw/1-s2.0-S1097276512008933-mmc2.xlsx"))
Hannum_CpGs <- Hannum_CpGs %>%
  select(Marker, Coefficient)
usethis::use_data(Hannum_CpGs, overwrite = TRUE)

#HRSInChPhenoAge Prep
HRSInCHPhenoAge_CpGs <-readr::read_csv("data-raw/HRSInChPhenoAge_CpG.csv")
HRSInCHPhenoAge_CpGs <- HRSInCHPhenoAge %>%
  select(CpG, Weight)
usethis::use_data(HRSInCHPhenoAge_CpGs, overwrite = TRUE)

#DunedinPACE and DunedinPoAm38 Prep
#   No prep needed, imported directly from github.com/danbelsky

#prcPhenoAge/ nonprcPhenoAge Prep
#   No prep needed, imported directly from github.com/MorganLevineLab

#prepare EpiToc, EpiToc2, HypoClock
load("data-raw/dataETOC2.Rd"); ## this loads the CpG information
hypoClock_CpGs <- dataETOC2.l[[3]]
EpiToc_CpGs <- dataETOC2.l[[2]]
EpiToc2_CpGs <- dataETOC2.l[[1]]
usethis::use_data(hypoClock_CpGs, overwrite = TRUE)
usethis::use_data(EpiToc_CpGs, overwrite = TRUE)
usethis::use_data(EpiToc2_CpGs, overwrite = TRUE)
rm(dataETOC2.l)

#prepare Garagnani
Garagnani_CpG <- "cg16867657"
usethis::use_data(Garagnani_CpG, overwrite = TRUE)

#prepare Bocklandt
Bocklandt_CpG <- "cg09809672"
usethis::use_data(Bocklandt_CpG, overwrite = TRUE)

#prepare Weidner
Weidner_CpGs <- c("cg02228185","cg25809905","cg17861230")
usethis::use_data(Weidner_CpGs, overwrite = TRUE)

#prepare Vidal-Bralo
VidalBralo_CpGs <- as.data.frame(readr::read_csv("data-raw/Clock-Vidal-Bralo-8-CpG.CSV"))
usethis::use_data(VidalBralo_CpGs, overwrite = TRUE)

#prepare Zhang
Zhang_10_CpG <-as.data.frame(readr::read_csv("data-raw/Clock-Zhang-10-CpG.CSV"))
usethis::use_data(Zhang_10_CpG, overwrite = TRUE)

#prepare Lin
Lin_CpGs <- as.data.frame(read.table("data-raw/Clock-Lin-99-CpG.txt",header=T))
usethis::use_data(Lin_CpGs, overwrite = TRUE)

#Prepare McCartney predictors
Smoking_CpGs <- as.data.frame(readr::read_csv("data-raw/DNAm_Smoking_McCartney.csv"))
BMI_CpGs <- as.data.frame(readr::read_csv("data-raw/DNAm_BMI_McCartney.csv"))
Alcohol_CpGs <- as.data.frame(readr::read_csv("data-raw/DNAm_Alcohol_McCartney.csv"))
usethis::use_data(Smoking_CpGs, overwrite = TRUE)
usethis::use_data(BMI_CpGs, overwrite = TRUE)
usethis::use_data(Alcohol_CpGs, overwrite = TRUE)

#prepare MiAge
MiAge_CpGs <- readr::read_csv("data-raw/MiAge_cpgs.csv")
MiAge_CpGs <- MiAge_CpGs %>%
  mutate(CpGs = CpG_site_ID) %>%
  select(CpGs, `Age-hyper/Age-hypo`)
usethis::use_data(MiAge_CpGs, overwrite = TRUE)
load("data-raw/MiAge_parameters.Rdata")
MiAge_parameters <- methyl.age
rm(methyl.age)
usethis::use_data(MiAge_parameters, overwrite = TRUE)

#prepare Zhang2019
Zhang2019_CpGs <- readr::read_delim("data-raw/en.coef",delim = " ")
colnames(Zhang2019_CpGs) <- c("CpG","coef")
#recall that intercept is 65.8
Zhang2019_CpGs <- Zhang2019_CpGs[-1,]
usethis::use_data(Zhang2019_CpGs, overwrite = TRUE)

#prepare PEDBE
PEDBE_CpGs <- readr::read_csv("data-raw/PEDBE_coef.csv")
PEDBE_CpGs <- PEDBE_CpGs %>%
  select(ID, Coef) %>%
  filter(ID != "(Intercept)")
usethis::use_data(PEDBE_CpGs, overwrite = TRUE)

#prepare DNAmTL
DNAmTL_CpGs <- readxl::read_xlsx(path = "~/Downloads/aging-v11i16-102173-supplementary-material-SD7.xlsx", skip = 5)
colnames(DNAmTL_CpGs)[1:2] <- c("ID","Coef")
DNAmTL_CpGs <- DNAmTL_CpGs %>%
  select(ID, Coef) %>%
  filter(ID != "Intercept")
usethis::use_data(DNAmTL_CpGs, overwrite = TRUE)

#prepare example Betas

load("data-raw/Ltd_CpG_beta_example.RData")
exampleBetas <- as.data.frame(exampleBetas)
exampleBetas <- exampleBetas[1:5,]
usethis::use_data(exampleBetas, overwrite = TRUE)

#prepare Horvath Online Calculator csv data
HorvathOnlineRef <- read.csv("data-raw/cgHorvathNew.csv")
HorvathOnlineRef <- HorvathOnlineRef %>%
  select(Name, Gene_ID, GenomeBuild, Chr, Accession)
usethis::use_data(HorvathOnlineRef, overwrite = TRUE)

#DNAmClockCortical
DNAmClockCortical_CpGs <- read.table("data-raw/CorticalClockCoefs.txt",stringsAsFactor=F,header=T)
colnames(DNAmClockCortical_CpGs) <- c("CpGs", "coef")
usethis::use_data(DNAmClockCortical_CpGs, overwrite = TRUE)

load("data-raw/Ref_DNAm_brain_values.rdat")
DNAmClockCortical_imputeRef <- ref
usethis::use_data(DNAmClockCortical_imputeRef, overwrite = TRUE)

#Knight gestational Age
Knight_CpGs <- readr::read_csv("data-raw/Knight_CpGs.csv")
colnames(Knight_CpGs) <- c("CpG","coef")
Knight_CpGs <- Knight_CpGs[-1,]
usethis::use_data(Knight_CpGs, overwrite = TRUE)

#Bohlin Gestational Age
load("data-raw/Bohlin_data.rda")
CpGs <- as.character(rownames(coef(UL.mod.cv, s = "lambda.min"))[as.logical(coef(UL.mod.cv, s = "lambda.min"))!=0][-1])
coef <- as.numeric(coef(UL.mod.cv, s = "lambda.min")[as.logical(coef(UL.mod.cv, s = "lambda.min"))!=0][-1])
Bohlin_CpGs <- data.frame(CpG = CpGs, coef = coef)
usethis::use_data(Bohlin_CpGs, overwrite = TRUE)

#Mayne Gestational Age
Mayne_CpGs <- readxl::read_xlsx("data-raw/Mayne_regression.xlsx", skip =1)
Mayne_CpGs <- Mayne_CpGs %>% select(Probe, Coefficient)
colnames(Mayne_CpGs) <- c("CpG","coef")
Mayne_CpGs <- Mayne_CpGs[-1,]
usethis::use_data(Mayne_CpGs, overwrite = TRUE)

Mayne_impute_table <- readxl::read_xlsx("data-raw/Mayne_meanTable.xlsx", skip = 1)
Mayne_impute <- as.vector(Mayne_impute_table$gold.mean)
names(Mayne_impute) <- Mayne_impute_table$probe
usethis::use_data(Mayne_impute, overwrite = TRUE)

#Lee placental clocks
Lee_all <- readr::read_csv("data-raw/Lee_raw_data.csv")
LeeRobust_CpGs <- data.frame(CpG = Lee_all$CpGs[Lee_all$Coefficient_RPC != 0][-1],
                             coef = Lee_all$Coefficient_RPC[Lee_all$Coefficient_RPC != 0][-1])
usethis::use_data(LeeRobust_CpGs, overwrite = TRUE)
LeeRefinedRobust_CpGs <- data.frame(CpG = Lee_all$CpGs[Lee_all$Coefficient_refined_RPC != 0][-1],
                             coef = Lee_all$Coefficient_refined_RPC[Lee_all$Coefficient_refined_RPC != 0][-1])
usethis::use_data(LeeRefinedRobust_CpGs, overwrite = TRUE)
LeeControl_CpGs <- data.frame(CpG = Lee_all$CpGs[Lee_all$Coefficient_CPC != 0][-1],
                             coef = Lee_all$Coefficient_CPC[Lee_all$Coefficient_CPC != 0][-1])
usethis::use_data(LeeControl_CpGs, overwrite = TRUE)









