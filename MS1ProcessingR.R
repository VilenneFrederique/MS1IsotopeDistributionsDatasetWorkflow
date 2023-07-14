MS1_Isotope_Distribution_Function <- function(Results){
  # Loading the required packages
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if(!require(BRAIN)){
    BiocManager::install("BRAIN")
    library(BRAIN)
  }
  if(!require(tidyverse)){
    install.packages("tidverse")
    library(tidyverse)
  }
  # BRAIN
  Results$NAs <- rowSums(is.na(Results))
  Results %>%
    mutate(Distribution = if_else(NAs >= 10, FALSE, TRUE )) -> Results
  Results["SpectralAngle"] <- NA
  Results["Carbons"] <- NA
  Results["Hydrogens"] <- NA
  Results["Oxygens"] <- NA
  Results["Nitrogens"] <- NA
  Results["Sulphurs"] <- NA
  Results["BRAINRelativeIsotopePeak1Intensity"] <- NA
  Results["BRAINRelativeIsotopePeak2Intensity"] <- NA
  Results["BRAINRelativeIsotopePeak3Intensity"] <- NA
  Results["BRAINRelativeIsotopePeak4Intensity"] <- NA
  Results["BRAINRelativeIsotopePeak5Intensity"] <- NA
  Results["BRAINRelativeIsotopePeak6Intensity"] <- NA
  for(index in 1:nrow(Results)){
    
    # Extract information from dataframe required
    Peptide <- Results[[index, "PeptideSequence"]]
    Modifications_Peptide <- Results[[index, "Modifications"]]
    CS <- Results[[index, "ChargeState"]]
    ExperimentalIsotopeDistr <- c(Results[[index, "IsotopePeak1Intensity"]],
                                  Results[[index, "IsotopePeak2Intensity"]],
                                  Results[[index, "IsotopePeak3Intensity"]],
                                  Results[[index, "IsotopePeak4Intensity"]],
                                  Results[[index, "IsotopePeak5Intensity"]],
                                  Results[[index, "IsotopePeak6Intensity"]])
    
    # BRAIN
    if(is.na(Modifications_Peptide)){
      Cysteines <- 0
      Methionines <- 0
    } else {
      Cysteines <- str_count(string = Modifications_Peptide , pattern = "(57.0214)")
      Methionines <- str_count(string = Modifications_Peptide , pattern = "(15.9949)")
    }
    AAComp <- getAtomsFromSeq(Peptide)
    ## Carbamidomethylation
    AAComp[["C"]] <- AAComp[["C"]] + 2 * Cysteines
    AAComp[["H"]] <- AAComp[["H"]] + 3 * Cysteines
    AAComp[["N"]] <- AAComp[["N"]] + 1 * Cysteines
    AAComp[["O"]] <- AAComp[["O"]] + 1 * Cysteines
    ## Oxidation
    AAComp[["O"]] <- AAComp[["O"]] + 1 * Methionines
    ## Hydrogens based on charge state
    AAComp[["H"]] <- AAComp[["H"]] + CS
    ## Brain
    BRAINResults <- useBRAIN(aC = AAComp, stopOption="nrPeaks", nrPeaks = 6)
    TheoreticalIsotopeDistr <- BRAINResults[["isoDistr"]]
    
    ## Append results to list
    Results[index, "Carbons"] <- AAComp[["C"]]
    Results[index, "Hydrogens"] <- AAComp[["H"]]
    Results[index, "Oxygens"] <- AAComp[["O"]]
    Results[index, "Nitrogens"] <- AAComp[["N"]]
    Results[index, "Sulphurs"] <- AAComp[["S"]]
    Results["BRAINRelativeIsotopePeak1Intensity"] <- BRAINResults[["isoDistr"]][1]
    Results["BRAINRelativeIsotopePeak2Intensity"] <- BRAINResults[["isoDistr"]][2]
    Results["BRAINRelativeIsotopePeak3Intensity"] <- BRAINResults[["isoDistr"]][3]
    Results["BRAINRelativeIsotopePeak4Intensity"] <- BRAINResults[["isoDistr"]][4]
    Results["BRAINRelativeIsotopePeak5Intensity"] <- BRAINResults[["isoDistr"]][5]
    Results["BRAINRelativeIsotopePeak6Intensity"] <- BRAINResults[["isoDistr"]][6]
    
    if(Results[[index, "Distribution"]] == TRUE){
      # Spectral Angle
      ## Remove NAs
      Non_NA_values <- which(!is.na(ExperimentalIsotopeDistr))
      ExperimentalIsotopeDistr <- ExperimentalIsotopeDistr[Non_NA_values]
      Theoretical_int <- TheoreticalIsotopeDistr[Non_NA_values]
      ## Normalize experimental data to 1
      Experimental_int <- ExperimentalIsotopeDistr / sum(ExperimentalIsotopeDistr, na.rm = TRUE)
      ## Compute SA
      Numerator <- sum(Theoretical_int * Experimental_int)
      Denominator <- sqrt(sum(Theoretical_int^2)) * sqrt(sum(Experimental_int^2))
      Spectral_angle <- acos(Numerator / Denominator)
      Results[index, "SpectralAngle"] <- Spectral_angle
    }
  }
  
  # Calculating the amount of detected peaks
  Results["NPeaks"] <- NA
  for(index in 1:nrow(Results)){
    Counter = 0
    if(!is.na(Results[[index, "IsotopePeak1Intensity"]])){
      Counter = Counter + 1
    }
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]])){
      Counter = Counter + 1
    }
    if(!is.na(Results[[index, "IsotopePeak3Intensity"]])){
      Counter = Counter + 1
    }
    if(!is.na(Results[[index, "IsotopePeak4Intensity"]])){
      Counter = Counter + 1
    }
    if(!is.na(Results[[index, "IsotopePeak5Intensity"]])){
      Counter = Counter + 1
    }
    if(!is.na(Results[[index, "IsotopePeak6Intensity"]])){
      Counter = Counter + 1
    }
    Results[index, "NPeaks"] <- Counter
  }
  
  # Looking if the peaks are consecutive or not
  Results["ConsecutivePeaks"] <- NA
  for(index in 1:nrow(Results)){
    NPeaks = Results[[index, "NPeaks"]]
    if (NPeaks == 1) {
      Results[index, "ConsecutivePeaks"] <- FALSE
    } else if(NPeaks == 2){
      if(!is.na(Results[[index, "IsotopePeak2Intensity"]])){
        Results[index, "ConsecutivePeaks"] <- TRUE
      } else {
        Results[index, "ConsecutivePeaks"] <- FALSE
      }
    }  else if(NPeaks == 3){
      if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]])){
        Results[index, "ConsecutivePeaks"] <- TRUE
      } else {
        Results[index, "ConsecutivePeaks"] <- FALSE
      }
    } else if(NPeaks == 4){
      if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]])){
        Results[index, "ConsecutivePeaks"] <- TRUE
      } else {
        Results[index, "ConsecutivePeaks"] <- FALSE
      }
    } else if(NPeaks == 5){
      if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]]) & !is.na(Results[[index, "IsotopePeak5Intensity"]])){
        Results[index, "ConsecutivePeaks"] <- TRUE
      } else {
        Results[index, "ConsecutivePeaks"] <- FALSE
      }
    } else if(NPeaks == 6){
      if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]]) & !is.na(Results[[index, "IsotopePeak5Intensity"]]) & !is.na(Results[[index, "IsotopePeak6Intensity"]])){
        Results[index, "ConsecutivePeaks"] <- TRUE
      } else {
        Results[index, "ConsecutivePeaks"] <- FALSE
      }
    }
  }
  
  Results %>%
    select(-NAs) -> Results
  
  return(Results)
}