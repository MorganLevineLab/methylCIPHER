
-   [MorganLevineLab/methylCIPHER](#morganlevinelabmethylcipher)
    -   [Installation](#installation)
    -   [Calculating Epigenetic Clocks and
        Predictors](#calculating-epigenetic-clocks-and-predictors)
        -  [Running single “clock” calculations] (#Running-single-“clock”-calculations)
        -  [Categories of Epigenetic Clocks] (#Categories-of-Epigenetic-clocks)
-   [Citations!](#citations)
    -   [citeMyClocks()](#citemyclocks)
    -   [calcCoreClocks Function](#calccoreclocks-function)
    -   [Chronological Age Predictors](#chronological-age-predictors-1)
    -   [Cancer and Mitotic Rates](#cancer-and-mitotic-rates-1)
    -   [Gestational & Pediatric Age](#gestational--pediatric-age)
    -   [Biological Age and Mortality
        Predictors](#biological-age-and-mortality-predictors-1)
    -   [Trait Predictors](#trait-predictors-1)
    -   [Bespoke Clocks](#bespoke-clocks)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# MorganLevineLab/methylCIPHER

<!-- badges: start -->
<!-- badges: end -->

The goal of methylCIPHER is to allow users to easily calculate their
choice of CpG clocks using simple commands, from a single source. CpG
epigenetic clocks are currently found in many places, and some require
users to send data to external portals, which is not advisable when
working with protected or restricted data. The current package allows
you to calculate reported epigenetic clocks, or where not precisely
disclosed, our best estimates–all performed locally on your own machine
in R Studio. We would like to acknowledge the authors of the original
clocks and their valuable contributions to aging research, and the
requisite citations for their clocks can be found at the bottom of the
current page. Please do not forget to cite them in your work!

## Installation

You can install the released version of methylCIPHER and its imported
packages from [Github](https://Cgithub.com/MorganLevineLab/methylCIPHER)
with:

``` r
devtools::install_github("MorganLevineLab/prcPhenoAge")
devtools::install_github("danbelsky/DunedinPoAm38")
devtools::install_github("MorganLevineLab/methylCIPHER")
```

## Calculating Epigenetic Clocks and Predictors

### Running single “clock” calculations

The current package contains a large number of currently available CpG
clocks or CpG based predictors. While we strove to be inclusive of such
published CpG-based epigenetic clocks to our knowledge, if you find we
are missing a clock, please contact us and we will do our best to
promptly include it, if possible. You can do so by raising an issue on
this repo or emailing us directly at
\<kyra\[dot\]thrush\[at\]yale\[dot\]edu\>.

In order to calculate a CpG clock, you simply need to use the
appropriate function, typically named “calc\[ClockNameHere\]”. For
example:

``` r
library(methylCIPHER)
calcPhenoAge(exampleBetas, examplePheno, imputation = F)
```

| name              | geo_accession | gender |  age | group | sample | PhenoAge |
|:------------------|:--------------|:-------|-----:|------:|-------:|---------:|
| 7786915023_R02C02 | GSM1343050    | M      | 57.9 |     1 |      1 | 52.29315 |
| 7786915135_R04C02 | GSM1343051    | M      | 42.0 |     1 |      2 | 41.05867 |
| 7471147149_R06C01 | GSM1343052    | M      | 47.4 |     1 |      3 | 43.54460 |
| 7786915035_R05C01 | GSM1343053    | M      | 49.3 |     1 |      4 | 43.96697 |
| 7786923035_R01C01 | GSM1343054    | M      | 52.5 |     1 |      5 | 40.35242 |

Alternatively, if you would just like to receive a vector with the clock
values to use, rather than appending it to an existing phenotype/
demographic dataframe, simply use:

``` r
calcPhenoAge(exampleBetas, imputation = F)
#> 7786915023_R02C02 7786915135_R04C02 7471147149_R06C01 7786915035_R05C01 
#>          52.29315          41.05867          43.54460          43.96697 
#> 7786923035_R01C01 
#>          40.35242
```

### Categories of Epigenetic Clocks

Due to the abundance of epigenetic clocks are overlapping reasons for
use, it is important to keep track of which clocks are most related to
each other. This will allow you to steer clear of multiple testing and
collinearity problems if you use all available clocks for your analysis.

Another important note is that with a few exceptions, the list of
following clocks were trained and validated almost exclusively in blood.
This can lead to a number of observed effects, which include
unreasonable shifts in age (intercept)As bespoke clocks are developed
for use in additional human tissues, we will include these in their own
section below.

#### Chronological Age Predictors

-   Bocklandt (calcBocklandt)
-   Garagnani (calcGaragnani)
-   **Hannum (calcHannum)**
-   **Horvath MultiTissue (calcHorvath1)**
-   Lin 99 CpG clock (calcLin)
-   Vidal-Bralo (calcVidalBralo)
-   Weidner (calcWeidner)
-   Zhang 10 CpG clock (calcZhang)

#### Cancer and Mitotic Rates

-   Epigenetic Timers of Cancer:
    -   EpiTOC (calcEpiTOC)
    -   **EpiTOC2 (calcEpiTOC2)**
-   solo-WCGWsites, or hypoclock (calcHypoClock)
-   MiAge Calculator (calcMiAge)

#### Gestational and Pediatric Age Predictors

-   Bohlin (calcBohlin)  
-   Knight (calcKnight)
-   Lee Gestational Age:
    -   calcLeeControl  
    -   calcLeeRobust  
    -   calcLeeRefinedRobust  
-   Mayne (calcMayne)  
-   Pediatric-Buccal-Epigenetic Clock (calcPEDBE)

#### Biological Age and Mortality Predictors

-   Dunedin Pace of Aging (DunedinPoAm38)
    -   Updated DunedinPACE: available
        [elsewhere](https://github.com/danbelsky/DunedinPACE)
-   **PhenoAge**:
    -   Original 513 CpG PhenoAge (calcPhenoAge)
        -   Polycomb Repressor Complex Associated Subscore
            ([calcprcPhenoAge](https://github.com/MorganLevineLab/prcPhenoAge))
        -   PRC dissociated Subscore (calcnonprcPhenoAge)
    -   Updated PhenoAge (calcHRSInChPhenoAge)
-   **GrimAge:** not available, use **formatHorvathOnline** [Horvath
    Online Calculator](http://dnamage.genetics.ucla.edu/home)

#### Trait Predictors

-   McCartney’s Predictors:
    -   Alcohol (calcAlcoholMcCartney)
    -   BMI (calcBMIMcCartney)
    -   Smoking (calcSmokingMcCartney)
-   PhosphatidylEthanol: Hazardous Alcohol Use (calcPEth)

#### Non-Blood Clocks

-   DNAmAge<sub>Cortical</sub>
    -   *A brain tissue trained clock.*
-   Horvath Skin & Blood (calcHorvath2)

### Running Multiple Epigenetic Clocks Simultaneously

#### Running the Core Epigenetic Clocks

There are a number of epigenetic clocks available in the literature, but
the majority of studies currently focus upon comparing a select few.
Therefore, we have provided a helper function to quickly calculate these
core clocks in just one line of code.

``` r
calcCoreClocks(exampleBetas, examplePheno)
#> Please remember to cite the core Clocks you have used! Please refer to the README.md file for assistance.
```

| name              | geo_accession | gender |  age | group | sample |   Hannum | Horvath1 | PhenoAge |  epiTOC2 |
|:------------------|:--------------|:-------|-----:|------:|-------:|---------:|---------:|---------:|---------:|
| 7786915023_R02C02 | GSM1343050    | M      | 57.9 |     1 |      1 | 62.33422 | 56.20544 | 52.29315 | 5012.412 |
| 7786915135_R04C02 | GSM1343051    | M      | 42.0 |     1 |      2 | 53.90139 | 47.45754 | 41.05867 | 4622.625 |
| 7471147149_R06C01 | GSM1343052    | M      | 47.4 |     1 |      3 | 53.59733 | 48.24172 | 43.54460 | 2956.300 |
| 7786915035_R05C01 | GSM1343053    | M      | 49.3 |     1 |      4 | 59.53767 | 51.53588 | 43.96697 | 3446.410 |
| 7786923035_R01C01 | GSM1343054    | M      | 52.5 |     1 |      5 | 58.80107 | 45.05448 | 40.35242 | 3245.157 |

#### Running A User-Defined List of Epigenetic Clocks

The user is welcome to specify a vector of clocks that they would like
to calculate, rather than running each individual clock calculation. In
this case, you will need to choose from the following options:

``` r
clockOptions()
#>  [1] "calcAlcoholMcCartney"            "calcBMIMcCartney"               
#>  [3] "calcBocklandt"                   "calcBohlin"                     
#>  [5] "calcDNAmClockCortical"           "calcDNAmTL"                     
#>  [7] "calcDunedinPoAm38"               "calcEpiTOC"                     
#>  [9] "calcEpiTOC2"                     "calcGaragnani"                  
#> [11] "calcHannum"                      "calcHorvath1"                   
#> [13] "calcHorvath2"                    "calcHRSInChPhenoAge"            
#> [15] "calcHypoClock"                   "calcKnight"                     
#> [17] "calcLeeControl"                  "calcLeeRefinedRobust"           
#> [19] "calcLeeRobust"                   "calcLin"                        
#> [21] "calcMayne"                       "calcMiAge"                      
#> [23] "calcPEDBE"                       "calcPhenoAge"                   
#> [25] "calcSmokingMcCartney"            "calcVidalBralo"                 
#> [27] "calcWeidner"                     "calcZhang"                      
#> [29] "calcZhang2019"                   "prcPhenoAge::calcPRCPhenoAge"   
#> [31] "prcPhenoAge::calcnonPRCPhenoAge" "DunedinPoAm38::PoAmProjector"
```

To do so, here is an example:

``` r
userClocks <- c("calcSmokingMcCartney","calcPhenoAge","calcEpiTOC2")
calcUserClocks(userClocks, exampleBetas, examplePheno, imputation = F)
```

| name              | geo_accession | gender |  age | group | sample | Smoking_McCartney | PhenoAge |  epiTOC2 |
|:------------------|:--------------|:-------|-----:|------:|-------:|------------------:|---------:|---------:|
| 7786915023_R02C02 | GSM1343050    | M      | 57.9 |     1 |      1 |          3.993508 | 52.29315 | 5012.412 |
| 7786915135_R04C02 | GSM1343051    | M      | 42.0 |     1 |      2 |          4.501657 | 41.05867 | 4622.625 |
| 7471147149_R06C01 | GSM1343052    | M      | 47.4 |     1 |      3 |          3.173744 | 43.54460 | 2956.300 |
| 7786915035_R05C01 | GSM1343053    | M      | 49.3 |     1 |      4 |          3.216788 | 43.96697 | 3446.410 |
| 7786923035_R01C01 | GSM1343054    | M      | 52.5 |     1 |      5 |          4.414541 | 40.35242 | 3245.157 |

### Missing beta values

Of course, all of the CpG clocks work best when you have all of the
necessary probes’ beta values for each sample. However, sometimes after
preprocessing, beta values will be removed for a variety of reasons. For
each CpG clock, you have the option to impute missing values for CpGs
that were removed across all samples. In this case, you will need to
impute using a vector of your choice (e.g. mean methylation values
across CpGs from an independent tissue-matched dataset). However, by
default, imputation will not be performed and the portion of the clock
that is reliant upon those CpGs will not be considered. To check quickly
whether this is the case for your data and clock(s) of interest, we have
created the following helper function:

``` r
getClockProbes(exampleBetas)
```

| Clock             | Total.Probes | Present.Probes | Percent.Present |
|:------------------|-------------:|---------------:|:----------------|
| Alcohol           |          450 |            450 | 100%            |
| BMI               |         1109 |           1109 | 100%            |
| Bocklandt         |            1 |              1 | 100%            |
| Bohlin            |          251 |              8 | 3%              |
| DNAmClockCortical |          347 |             33 | 10%             |
| DNAmTL            |          140 |             11 | 8%              |
| EpiToc2           |          163 |            163 | 100%            |
| EpiToc            |          385 |            385 | 100%            |
| Garagnani         |            1 |              1 | 100%            |
| HRSInCHPhenoAge   |          959 |            959 | 100%            |
| Hannum            |           71 |             71 | 100%            |
| Horvath1          |          353 |            353 | 100%            |
| Horvath2          |          391 |            390 | 100%            |
| Knight            |            1 |              0 | 0%              |
| LeeControl        |          546 |             13 | 2%              |
| LeeRefinedRobust  |          395 |              9 | 2%              |
| LeeRobust         |          558 |              9 | 2%              |
| Lin               |           99 |             39 | 39%             |
| Mayne             |            1 |              0 | 0%              |
| MiAge             |          268 |              4 | 1%              |
| PEDBE             |           94 |             14 | 15%             |
| PhenoAge          |          513 |            513 | 100%            |
| Smoking           |          233 |            233 | 100%            |
| VidalBralo        |            8 |              5 | 62%             |
| Weidner           |            3 |              3 | 100%            |
| Zhang2019         |          514 |            131 | 25%             |
| Zhang             |           10 |             10 | 100%            |
| hypoClock         |          678 |            678 | 100%            |

Please note that this will not count columns of NAs for named CpGs as
missing! If you want to check for this you can run the following line of
code to find the column numbers that are all NAs. If you get “named
integer(0)” then you don’t have any. We recommend that you remove any
identified columns from your beta matrix entirely to avoid errors, and
then rerun the code producing the table above.

``` r
which(apply(exampleBetas, 2, function(x)all(is.na(x))))
```

In the case that you have CpGs missing from only some samples, we
encourage you to be aware of this early on. Run the following line and
check that it is 0.

``` r
sum(is.na(betaMatrix))
```

If this does not end up being 0, you might consider running mean
imputation within your data so that NA values for single/ few samples at
least have mean values rather than being ignored.

### Alternate Methods

#### PC Clocks

If you would like to lower the noise in your data, which is helpful for
applications in clinical, intervention, or longitudinal data in
particular, you might consider instead using our recently developed PC
Clocks. You can use [this
code](https://github.com/MorganLevineLab/PC-Clocks) and refer to [this
publication](insert-paper-doi-when-published).

#### Modules Clocks

\[should we include this?\]

# Citations!

The current package is intended to increase visibility and accesibility
for the many epigenetic clocks currently in circulation. Please do not
forget to cite not only this work, but the works you use in your
research the correspond to the appropriate epigenetic clocks.

## citeMyClocks()

As an alternative to manually checking the list below, you can simply
use the function **citeMyClocks** to automatically print out APA style
citations of the appropriate publications. You can use a character
vector or input with the function name(s) you used, or the same
*userClocks* list you may have used for **calcUserClocks**.

``` r
citeMyClocks(userClocks, prettyprint = TRUE)
#> Clock Citations:
#> Levine, M. E., Lu, A. T., Quach, A., Chen, B. H., Assimes, T. L., Bandinelli, S., … Horvath, S. (2018). 
#>  An epigenetic biomarker of aging for lifespan and healthspan. 
#>  Aging, 10(4), 573–591. https://doi.org/10.18632/aging.101414
#> McCartney, D. L., Hillary, R. F., Stevenson, A. J., Ritchie, S. J., Walker, R. M., Zhang, Q., … Marioni, R. E. (2018). 
#>  Epigenetic prediction of complex traits and death. 
#>  Genome Biology, 19(1), 136. https://doi.org/10.1186/s13059-018-1514-1
#> Teschendorff, A. E. (2020). 
#>  A comparison of epigenetic mitotic-like clocks for cancer risk prediction. 
#>  Genome Medicine, 12(1), 56. https://doi.org/10.1186/s13073-020-00752-3
```

*Please note that by default the function formats the output for nicer
console output. Change to prettyprint = FALSE for a vector you can save
in your workspace.*

## calcCoreClocks Function

If you used the **calcCoreClocks** function then the list of citations
you will need is as below. Otherwise, please find the appropriate
citations for the calculations you did run in their respective sections
below.

-   Hannum, G., Guinney, J., Zhao, L., Zhang, L., Hughes, G., Sadda, S.
    V., … Zhang, K. (2013). Genome-wide Methylation Profiles Reveal
    Quantitative Views of Human Aging Rates. Molecular Cell, 49(2),
    359–367. <https://doi.org/10.1016/j.molcel.2012.10.016>  
-   Horvath, S. (2013). DNA methylation age of human tissues and cell
    types. Genome Biology, 14(10), 3156.
    <https://doi.org/10.1186/gb-2013-14-10-r115>  
-   Levine, M. E., Lu, A. T., Quach, A., Chen, B. H., Assimes, T. L.,
    Bandinelli, S., … Horvath, S. (2018). An epigenetic biomarker of
    aging for lifespan and healthspan. Aging, 10(4), 573–591.
    <https://doi.org/10.18632/aging.101414>  
-   Teschendorff, A. E. (2020). A comparison of epigenetic mitotic-like
    clocks for cancer risk prediction. Genome Medicine, 12(1), 56.
    <https://doi.org/10.1186/s13073-020-00752-3>

## Chronological Age Predictors

-   **calcBocklandt:** Bocklandt, S., Lin, W., Sehl, M. E., Sánchez, F.
    J., Sinsheimer, J. S., Horvath, S., & Vilain, E. (2011). Epigenetic
    Predictor of Age. PLOS ONE, 6(6), e14821.
    <https://doi.org/10.1371/journal.pone.0014821>  
-   **calcGaragnani:** Garagnani, P., Bacalini, M. G., Pirazzini, C.,
    Gori, D., Giuliani, C., Mari, D., … Franceschi, C. (2012).
    Methylation of ELOVL2 gene as a new epigenetic marker of age. Aging
    Cell, 11(6), 1132–1134. <https://doi.org/10.1111/acel.12005>  
-   **calcHannum:** Hannum, G., Guinney, J., Zhao, L., Zhang, L.,
    Hughes, G., Sadda, S. V., … Zhang, K. (2013). Genome-wide
    Methylation Profiles Reveal Quantitative Views of Human Aging Rates.
    Molecular Cell, 49(2), 359–367.
    <https://doi.org/10.1016/j.molcel.2012.10.016>  
-   **calcHorvath1:** Horvath, S. (2013). DNA methylation age of human
    tissues and cell types. Genome Biology, 14(10), 3156.
    <https://doi.org/10.1186/gb-2013-14-10-r115>  
-   **calcLin:** Lin, Q., Weidner, C. I., Costa, I. G., Marioni, R. E.,
    Ferreira, M. R. P., Deary, I. J., & Wagner, W. (2016). DNA
    methylation levels at individual age-associated CpG sites can be
    indicative for life expectancy. Aging, 8(2), 394–401.
    <https://doi.org/10.18632/aging.100908>  
-   **calcVidalBralo:** Vidal-Bralo, L., Lopez-Golan, Y., & Gonzalez, A.
    (2016). Simplified assay for epigenetic age estimation in whole
    blood of adults. Frontiers in Genetics, 7(JUL), 1–7.
    <https://doi.org/10.3389/fgene.2016.00126>  
-   **calcWeidner:** Weidner, C. I., Lin, Q., Koch, C. M., Eisele, L.,
    Beier, F., Ziegler, P., … Wagner, W. (2014). Aging of blood can be
    tracked by DNA methylation changes at just three CpG sites. Genome
    Biology, 15(2). <https://doi.org/10.1186/gb-2014-15-2-r24>  
-   **calcZhang:** Zhang, Y., Wilson, R., Heiss, J., Breitling, L. P.,
    Saum, K. U., Schöttker, B., … Brenner, H. (2017). DNA methylation
    signatures in peripheral blood strongly predict all-cause mortality.
    Nature Communications, 8. <https://doi.org/10.1038/ncomms14617>  
-   **calcZhang2019:** Zhang, Q., Vallerga, C. L., Walker, R. M., Lin,
    T., Henders, A. K., Montgomery, G. W., … Visscher, P. M. (2019).
    Improved precision of epigenetic clock estimates across tissues and
    its implication for biological ageing. Genome Medicine, 11(1), 54.
    <https://doi.org/10.1186/s13073-019-0667-1>

## Cancer and Mitotic Rates

-   **calcDNAmTL:** Lu, A. T., Seeboth, A., Tsai, P. C., Sun, D., Quach,
    A., Reiner, A. P., … Horvath, S. (2019). DNA methylation-based
    estimator of telomere length. Aging, 11(16), 5895–5923.
    <https://doi.org/10.18632/aging.102173>  
-   **calcEpiTOC:** Yang, Z., Wong, A., Kuh, D., Paul, D. S., Rakyan, V.
    K., Leslie, R. D., … Teschendorff, A. E. (2016). Correlation of an
    epigenetic mitotic clock with cancer risk. Genome Biology, 17(1),
    1–18. <https://doi.org/10.1186/s13059-016-1064-3>  
-   **calcEpiTOC2 + calcHypoClock:** Teschendorff, A. E. (2020). A
    comparison of epigenetic mitotic-like clocks for cancer risk
    prediction. Genome Medicine, 12(1), 56.
    <https://doi.org/10.1186/s13073-020-00752-3>  
-   **calcMiAge:** Youn, A., & Wang, S. (2018). The MiAge Calculator: a
    DNA methylation-based mitotic age calculator of human tissue types.
    Epigenetics, 13(2), 192–206.
    <https://doi.org/10.1080/15592294.2017.1389361>

## Gestational & Pediatric Age

-   **calcBohlin:** Bohlin, J., Håberg, S. E., Magnus, P., Reese, S. E.,
    Gjessing, H. K., Magnus, M. C., … Nystad, W. (2016). Prediction of
    gestational age based on genome-wide differentially methylated
    regions. Genome Biology, 17(1), 1–9.
    <https://doi.org/10.1186/s13059-016-1063-4>  
-   **calcKnight:** Knight, A. K., Craig, J. M., Theda, C.,
    Bækvad-Hansen, M., Bybjerg-Grauholm, J., Hansen, C. S., …
    Smith, A. K. (2016). An epigenetic clock for gestational age at
    birth based on blood methylation data. Genome Biology, 17(1), 1–11.
    <https://doi.org/10.1186/s13059-016-1068-z>  
-   Lee, Y., Choufani, S., Weksberg, R., Wilson, S. L., Yuan, V., Burt,
    A., … Horvath, S. (2019). Placental epigenetic clocks: Estimating
    gestational age using placental DNA methylation levels. Aging,
    11(12), 4238–4253. <https://doi.org/10.18632/aging.102049>
    -   **calcLeeControl:** Control Placental Sample Age Calculator  
    -   **calcLeeRobust:** Robust predictor of Placental Age  
    -   **calcLeeRefinedRobust:** A refined version of the Robust
        Predictor of Placental Age  
-   **calcMayne:** Mayne, B. T., Leemaqz, S. Y., Smith, A. K., Breen,
    J., Roberts, C. T., & Bianco-Miotto, T. (2017). Accelerated
    placental aging in early onset preeclampsia pregnancies identified
    by DNA methylation. Epigenomics, 9(3), 279–289.
    <https://doi.org/10.2217/epi-2016-0103>  
-   **calcPEDBE:** McEwen, L. M., O’Donnell, K. J., McGill, M. G.,
    Edgar, R. D., Jones, M. J., MacIsaac, J. L., … Kobor, M. S. (2020).
    The PedBE clock accurately estimates DNA methylation age in
    pediatric buccal cells. PNAS, 117(38), 23329–23335.
    <https://doi.org/10.1073/pnas.1820843116>

## Biological Age and Mortality Predictors

-   **calcDunedinPoAm38:** Belsky, D. W., Caspi, A., Arseneault, L.,
    Baccarelli, A., Corcoran, D. L., Gao, X., … Moffitt, T. E. (2020).
    Quantification of the pace of biological aging in humans through a
    blood test, the DunedinPoAm DNA methylation algorithm. ELife, 9,
    e54870. <https://doi.org/10.7554/eLife.54870>  
-   **dunedinPACE (DunedinPACE::PoAmProjector())** Belsky, D. W., Caspi,
    A., Corcoran, D. L., Sugden, K., Poulton, R., Arseneault, L., …
    Moffitt, T. E. (2022). DunedinPACE, a DNA methylation biomarker of
    the pace of aging. ELife, 11, e73420.
    <https://doi.org/10.7554/eLife.73420>
-   **calcPhenoAge:** Levine, M. E., Lu, A. T., Quach, A., Chen, B. H.,
    Assimes, T. L., Bandinelli, S., … Horvath, S. (2018). An epigenetic
    biomarker of aging for lifespan and healthspan. Aging, 10(4),
    573–591. <https://doi.org/10.18632/aging.101414>  
-   **prcPhenoAge::calcprcPhenoAge() / calcnonprcPhenoAge()** Add
    prcPhenoAge manuscript citation
-   **calcHRSInChPhenoAge:** \[cite PC Clocks paper here\]
-   *GrimAge:* Lu, A. T., Quach, A., Wilson, J. G., Reiner, A. P., Aviv,
    A., Raj, K., … Horvath, S. (2019). DNA methylation GrimAge strongly
    predicts lifespan and healthspan. Aging, 11(2), 303–327.
    <https://doi.org/10.18632/aging.101684>

## Trait Predictors

-   McCartney, D. L., Hillary, R. F., Stevenson, A. J., Ritchie, S. J.,
    Walker, R. M., Zhang, Q., … Marioni, R. E. (2018). Epigenetic
    prediction of complex traits and death. Genome Biology, 19(1), 136.
    <https://doi.org/10.1186/s13059-018-1514-1>
    -   **calcAlcoholMcCartney**
    -   **calcBMIMcCartney**
    -   **calcSmokingMcCartney**
-   **calcPEth:** Liang, X., Justice, A. C., So-Armah, K., Krystal, J.
    H., Sinha, R., & Xu, K. (2021). DNA methylation signature on
    phosphatidylethanol, not on self-reported alcohol consumption,
    predicts hazardous alcohol consumption in two distinct populations.
    Molecular Psychiatry, 26(6), 2238–2253.
    <https://doi.org/10.1038/s41380-020-0668-x>

## Bespoke Clocks

-   **calcDNAmClockCortical:** Shireby, G. L., Davies, J. P.,
    Francis, P. T., Burrage, J., Walker, E. M., Neilson, G. W. A., …
    Mill, J. (2020). Recalibrating the epigenetic clock: implications
    for assessing biological age in the human cortex. Brain, 1–13.
    <https://doi.org/10.1093/brain/awaa334>  
-   **calcPEDBE:** McEwen, L. M., O’Donnell, K. J., McGill, M. G.,
    Edgar, R. D., Jones, M. J., MacIsaac, J. L., … Kobor, M. S. (2020).
    The PedBE clock accurately estimates DNA methylation age in
    pediatric buccal cells. Proceedings of the National Academy of
    Sciences of the United States of America, 117(38), 23329–23335.
    <https://doi.org/10.1073/pnas.1820843116>  
-   **calcHorvath2:** Horvath, S., Oshima, J., Martin, G. M., Lu, A. T.,
    Quach, A., Cohen, H., … Raj, K. (2018). Epigenetic clock for skin
    and blood cells applied to Hutchinson Gilford Progeria Syndrome and
    ex vivo studies. Aging, 10(7), 1758–1775.
    <https://doi.org/10.18632/aging.101508>
