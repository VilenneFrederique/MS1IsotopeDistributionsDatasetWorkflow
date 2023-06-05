# Libraries
library(tidyverse)
library(BRAIN)
library(ggvenn)
library(readxl)
library(writexl)
library(cowplot)
library(caret)
library(xlsx)


# MSFragger results
PSMs_12042 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12042/psm.tsv")
Peptides_12042 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12042/peptide.tsv")
Protein_12042 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12042/protein.tsv")
PSMs_12043 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12043/psm.tsv")
Peptides_12043 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12043/peptide.tsv")
Protein_12043 <- read.delim("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/A11-12043/protein.tsv")

# Save MSFragger results to single Excel
write.xlsx(PSMs_12042, file="MSFraggerResultsUPS.xlsx", sheetName="PSMs_A11-12042", row.names=FALSE)
write.xlsx(Peptides_12042, file="MSFraggerResultsUPS.xlsx", sheetName="Peptides_A11-12042", append=TRUE, row.names=FALSE)
write.xlsx(Protein_12042, file="MSFraggerResultsUPS.xlsx", sheetName="Proteins_A11-12042", append=TRUE, row.names=FALSE)
write.xlsx(PSMs_12043, file="MSFraggerResultsUPS.xlsx", sheetName="PSMs_A11-12043", append=TRUE, row.names=FALSE)
write.xlsx(Peptides_12042, file="MSFraggerResultsUPS.xlsx", sheetName="Peptides_A11-12043", append=TRUE, row.names=FALSE)
write.xlsx(Protein_12042, file="MSFraggerResultsUPS.xlsx", sheetName="Proteins_A11-12043", append=TRUE, row.names=FALSE)

# Venn Diagrams
## Protein level
Protein_list_12042 <- Protein_12042$Protein.ID
Protein_list_12043 <- Protein_12043$Protein.ID
Proteins <-list("A11-12042"=Protein_list_12042, "A11-12043"=Protein_list_12043)
PlotVennProt <- ggvenn(data = Proteins,
                       fill_color = c("mediumblue", "red2"),
                       fill_alpha = 0.5,
                       text_size = 6)
PlotVennProt
Unique_Prot_12042 = setdiff(Protein_list_12042, Protein_list_12043)
Unique_Prot_12042
Unique_Prot_12043 = setdiff(Protein_list_12043, Protein_list_12042)
Unique_Prot_12043

## Peptide level
Peptide_list_12042 <- Peptides_12042$Peptide
Peptide_list_12043 <- Peptides_12043$Peptide
Peptides <-list("A11-12042"= Peptide_list_12042, "A11-12043"= Peptide_list_12043)
PlotVennPept <- ggvenn(data = Peptides,
                       fill_color = c("mediumblue", "red2"),
                       fill_alpha = 0.5,
                       text_size = 6)
PlotVennPept

## Looking into the uniquely identified peptides by both A11-12042 and A11-12043
Unique_Pept_12042 = setdiff(Peptide_list_12042, Peptide_list_12043)
Unique_Pept_12042
Unique_Pept_12043 = setdiff(Peptide_list_12043, Peptide_list_12042)
Unique_Pept_12043

# Percentage proteins detected
## Creating a dataframe with all proteins, their concentration and if they were detected in the sample or not
Proteins_matrix = matrix(nrow = 48, ncol = 4)
Proteins_df = data.frame(Proteins_matrix)
colnames(Proteins_df) <- c("Protein.ID", "Concentration", "A11-12042", "A11-12043")
All_proteins_list <- c("P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human)",
                       "Q15843ups|NEDD8_HUMAN_UPS NEDD8 (Chain 1-81) - Homo sapiens (Human)",
                       "P01112ups|RASH_HUMAN_UPS GTPase HRas (Chain 1-189) - Homo sapiens (Human)",
                       "P04040ups|CATA_HUMAN_UPS Catalase (Chain 2-527) - Homo sapiens (Human)",
                       "P02753ups|RETBP_HUMAN_UPS Retinol-binding protein 4 (Chain 19-201) - Homo sapiens (Human)",
                       "P06396ups|GELS_HUMAN_UPS Gelsolin (Chain 28-782) - Homo sapiens (Human)",
                       "P01375ups|TNFA_HUMAN_UPS Tumor necrosis factor, soluble form (Chain 77-233) - Homo sapiens (Human)",
                       "P15559ups|NQO1_HUMAN_UPS NAD(P)H dehydrogenase [quinone] 1 (Chain 2-274) - Homo sapiens (Human)",
                       "Q06830ups|PRDX1_HUMAN_UPS Peroxiredoxin 1 (Chain 2-199) - Homo sapiens (Human)",
                       "P00167ups|CYB5_HUMAN_UPS Cytochrome b5 (Chain 1-134, N-terminal His tag) - Homo sapiens (Human)",
                       "P06732ups|KCRM_HUMAN_UPS Creatine kinase M-type (Chain 1-381) - Homo sapiens (Human)",
                       "P02741ups|CRP_HUMAN_UPS C-reactive protein (Chain 19-224) - Homo sapiens (Human)",
                       "P61626ups|LYSC_HUMAN_UPS Lysozyme C (Chain 19-148) - Homo sapiens (Human)",
                       "P16083ups|NQO2_HUMAN_UPS Ribosyldihydronicotinamide dehydrogenase [quinone] (Chain 2-231) - Homo sapiens (Human)",
                       "P10145ups|IL8_HUMAN_UPS Interleukin-8, IL-8 (Chain 28-99) - Homo sapiens (Human)",
                       "P61769ups|B2MG_HUMAN_UPS Beta-2-microglobulin (Chain 21-119) - Homo sapiens (Human)",
                       "P02144ups|MYG_HUMAN_UPS Myoglobin (Chain 2-154) - Homo sapiens (Human)",
                       "P08263ups|GSTA1_HUMAN_UPS Glutathione S-transferase A1 (Chain 2-222) - Homo sapiens (Human)",
                       "P55957ups|BID_HUMAN_UPS BH3-interacting domain death agonist (Chain 1-195) - Homo sapiens (Human)",
                       "P00709ups|LALBA_HUMAN_UPS Alpha-lactalbumin (Chain 20-142) - Homo sapiens (Human)",
                       "P69905ups|HBA_HUMAN_UPS Hemoglobin subunit alpha (Chain 2-142) - Homo sapiens (Human)",
                       "P01344ups|IGF2_HUMAN_UPS Insulin-like growth factor II (Chain 25-91) - Homo sapiens (Human)",
                       "P00915ups|CAH1_HUMAN_UPS Carbonic anhydrase 1 (Chain 2-261) - Homo sapiens (Human)",
                       "P08758ups|ANXA5_HUMAN_UPS Annexin A5 (Chain 2-320) - Homo sapiens (Human)",
                       "P00918ups|CAH2_HUMAN_UPS Carbonic anhydrase 2 (Chain 2-260) - Homo sapiens (Human)",
                       "P62937ups|PPIA_HUMAN_UPS Peptidyl-prolyl cis-trans isomerase A (Chain 1-165, N terminal His tag)- Homo sapiens (Human)",
                       "P00441ups|SODC_HUMAN_UPS Superoxide dismutase [Cu-Zn] (Chain 2-154) - Homo sapiens (Human)",
                       "P05413ups|FABPH_HUMAN_UPS Fatty acid-binding protein, heart (Chain 2-133) - Homo sapiens (Human)",
                       "O00762ups|UBE2C_HUMAN_UPS Ubiquitin-conjugating enzyme E2 C (Chain 1-179, N-terminal His tag)- Homo sapiens (Human)",
                       "P41159ups|LEP_HUMAN_UPS Leptin (Chain 22-167) - Homo sapiens (Human)",
                       "P02788ups|TRFL_HUMAN_UPS Lactotransferrin (Chain 20-710) - Homo sapiens (Human)",
                       "P09211ups|GSTP1_HUMAN_UPS Glutathione S-transferase P (Chain 2-210) - Homo sapiens (Human)",
                       "P01031ups|CO5_HUMAN_UPS Complement C5 (C5a anaphylatoxin) (Chain 678-751) - Homo sapiens (Human)",
                       "P10636-8ups|TAU_HUMAN_UPS Microtubule-associated protein tau {Isoform Tau-F (Tau-4)} (Chain 2-441) - Homo sapiens (Human)",
                       "P01133ups|EGF_HUMAN_UPS Pro-Epidermal growth factor (EGF) (Chain 971-1023) - Homo sapiens (Human)",
                       "P63165ups|SUMO1_HUMAN_UPS Small ubiquitin-related modifier 1 (Chain 1-97, N-terminal GST tag) - Homo sapiens (Human)",
                       "P51965ups|UB2E1_HUMAN_UPS Ubiquitin-conjugating enzyme E2 E1 (Chain 1-193, N terminal His tag)- Homo sapiens (Human)",
                       "P12081ups|SYHC_HUMAN_UPS Histidyl-tRNA synthetase, cytoplasmic (Chain 1-509, C terminal His tag) - Homo sapiens (Human)",
                       "P99999ups|CYC_HUMAN_UPS Cytochrome c (Chain 2-105) - Homo sapiens (Human)",
                       "P02787ups|TRFE_HUMAN_UPS Serotransferrin (Chain 20-698) - Homo sapiens (Human)",
                       "P10599ups|THIO_HUMAN_UPS Thioredoxin (Chain 2-105, N-terminal His tag)- Homo sapiens (Human)",
                       "P62988ups|UBIQ_HUMAN_UPS Ubiquitin (Chain 1-76, N-terminal His tag) - Homo sapiens (Human)",
                       "P01008ups|ANT3_HUMAN_UPS Antithrombin-III (Chain 33-464) - Homo sapiens (Human)",
                       "P01127ups|PDGFB_HUMAN_UPS Platelet-derived growth factor B chain (Chain 82-190) - Homo sapiens (Human)",
                       "P68871ups|HBB_HUMAN_UPS Hemoglobin subunit beta (Chain 2-147) - Homo sapiens (Human)",
                       "P63279ups|UBC9_HUMAN_UPS SUMO-conjugating enzyme UBC9 (Chain 1-158) - Homo sapiens (Human)",
                       "O76070ups|SYUG_HUMAN_UPS Gamma-synuclein (Chain 1-127) - Homo sapiens (Human)",
                       "P01579ups|IFNG_HUMAN_UPS Interferon Gamma (Chain 23-166) - Homo sapiens (Human)")
All_concentrations <- c(50,
                        0.5,
                        0.005,
                        5,
                        0.5,
                        0.005,
                        0.0005,
                        5,
                        5,
                        5,
                        0.5,
                        0.0005,
                        0.5,
                        0.5,
                        0.0005,
                        0.05,
                        5,
                        0.05,
                        0.05,
                        0.05,
                        50,
                        0.05,
                        50,
                        0.0005,
                        50,
                        5,
                        0.0005,
                        0.0005,
                        0.005,
                        50,
                        0.0005,
                        0.05,
                        50,
                        0.0005,
                        5,
                        5,
                        0.005,
                        0.5,
                        0.005,
                        0.005,
                        0.05,
                        50,
                        0.05,
                        0.05,
                        50,
                        0.5,
                        0.05,
                        0.005)
Proteins_df$Protein.ID <- All_proteins_list
Proteins_df$Concentration <- All_concentrations
Proteins_df$`A11-12042` <- ifelse(Proteins_df$Protein.ID %in% Protein_list_12042, 1, 0)
Proteins_df$`A11-12043` <- ifelse(Proteins_df$Protein.ID %in% Protein_list_12043, 1, 0)

## Table of detected proteins by each sample
table(Proteins_df$Concentration, Proteins_df$`A11-12042`)
table(Proteins_df$Concentration, Proteins_df$`A11-12043`)

## Preparing to plot the percentages of Proteins detected in function of concentration per sample
Percentages_DF <- matrix(nrow = 6, ncol = 3)
Percentages_DF <- data.frame(Percentages_DF)
colnames(Percentages_DF) <- c("Concentrations", "A11-12042", "A11-12043")
Percentages_DF$Concentrations <- c(50, 5, 0.5, 0.05, 0.005, 0.0005)
Percentages_DF$`A11-12042` <- c(100, 100, 100, 90, 57.14286, 0)
Percentages_DF$`A11-12043` <- c(100, 100, 100, 90, 57.14286, 0)
Percentages_DF %>%
  gather(key = "Sample", value = "Percentages", 2:3) -> Percentages_DF
Percentages_DF %>%
  ggplot(aes(x = factor(Concentrations), y = Percentages, fill = Sample, colour = Sample)) +
  geom_col(alpha = 0.5, position = "dodge") +
  theme_minimal() +
  scale_fill_manual(values = c("mediumblue", "red2")) +
  scale_color_manual(values = c("mediumblue", "red2")) +
  xlab("Concentration (pmol)") +
  ylab("Detected (%)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) -> PlotPercProt
PlotPercProt
 
## Plotting the coverage proteins in function of concentration per sample
Protein_12042 %>%
  select(Protein.ID, Coverage) -> Protein_12042_subset
Proteins_df <- left_join(x = Proteins_df, y = Protein_12042_subset, by = "Protein.ID")
Proteins_df %>%
  dplyr::rename("Coverage_12042" = "Coverage") -> Proteins_df
Protein_12043 %>%
  select(Protein.ID, Coverage) -> Protein_12043_subset
Proteins_df <- left_join(x = Proteins_df, y = Protein_12043_subset, by = "Protein.ID")
Proteins_df %>%
  dplyr::rename("Coverage_12043" = "Coverage") -> Proteins_df
Proteins_df %>%
  select(Concentration, Coverage_12042, Coverage_12043) -> Coverages
Coverages %>%
  gather(key = "Sample", value = "Coverage", 2:3) -> Coverages
Coverages$Sample <- str_replace(Coverages$Sample, "Coverage_", "A11-")
Coverages <- Coverages %>% replace(is.na(.), 0)
Coverages %>%
  group_by(Concentration, Sample) %>%
  summarise(Average_Coverage = mean(Coverage),
            SD_Coverage = sd(Coverage)) -> Coverages2
Coverages %>% 
  ggplot(aes(x = factor(Concentration), y = Coverage, fill = Sample, colour = Sample)) +
  geom_boxplot(alpha = 0.5, position = "dodge") +
  #geom_jitter(alpha = 0.5, width = 0.5) +
  theme_minimal() +
  scale_fill_manual(values = c("mediumblue", "red2")) +
  scale_colour_manual(values = c("mediumblue", "red2")) +
  xlab("Concentration (pmol)") +
  ylab("Coverage (%)") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) -> PlotCovProt
PlotCovProt

## Plotting all plots together
ggdraw() +
  draw_plot(PlotVennPept, x = 0, y = .5, width = .5, height = 0.5) +
  draw_plot(PlotVennProt, x = 0.5, y = .5, width = .5, height = 0.5) +
  draw_plot(PlotPercProt, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(PlotCovProt, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.5, 0, 0.5), y = c(0.95, 0.95, 0.55, 0.55))


#####################################################################################
# Merging of Excel files
A11_12042_mzML <- read_excel("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/Article/A11-12042.mzMLObsMZ30Sec.xlsx")
A11_12043_mzML <- read_excel("G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/Article/A11-12043.mzMLObsMZ30Sec.xlsx")
Results <- bind_rows(A11_12042_mzML, A11_12043_mzML)


###################################################################################
# Computing the spectral angle as a quality measure
## BRAIN
Results$NAs <- rowSums(is.na(Results))
Results %>%
  mutate(Distribution = if_else(NAs >= 10, "No", "Yes" )) -> Results
Results["SpectralAngle"] <- NA
Results["Carbons"] <- NA
Results["Hydrogens"] <- NA
Results["Oxigens"] <- NA
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
  ## Brain
  BRAINResults <- useBRAIN(aC = AAComp, stopOption="nrPeaks", nrPeaks = 6)
  TheoreticalIsotopeDistr <- BRAINResults[["isoDistr"]]
  
  ## Append results to list
  Results[index, "Carbons"] <- AAComp[["C"]]
  Results[index, "Hydrogens"] <- AAComp[["H"]]
  Results[index, "Oxigens"] <- AAComp[["O"]]
  Results[index, "Nitrogens"] <- AAComp[["N"]]
  Results[index, "Sulphurs"] <- AAComp[["S"]]
  Results["BRAINRelativeIsotopePeak1Intensity"] <- BRAINResults[["isoDistr"]][1]
  Results["BRAINRelativeIsotopePeak2Intensity"] <- BRAINResults[["isoDistr"]][2]
  Results["BRAINRelativeIsotopePeak3Intensity"] <- BRAINResults[["isoDistr"]][3]
  Results["BRAINRelativeIsotopePeak4Intensity"] <- BRAINResults[["isoDistr"]][4]
  Results["BRAINRelativeIsotopePeak5Intensity"] <- BRAINResults[["isoDistr"]][5]
  Results["BRAINRelativeIsotopePeak6Intensity"] <- BRAINResults[["isoDistr"]][6]
  
  if(Results[[index, "Distribution"]] == "Yes"){
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

## Plotting the distribution of the spectral angle per sample
Results %>%
  mutate(RawFile = recode(RawFile, 
                          "A11-12042.mzML" = "A11-12042",
                          "A11-12043.mzML" = "A11-12043")) %>%
  ggplot(aes(x = SpectralAngle, fill = RawFile)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  scale_fill_manual(name = "Sample", values = c("mediumblue", "red2")) +
  facet_wrap(.~RawFile, ncol = 1) +
  xlab(label = "Spectral angle") +
  geom_vline(data = data.frame(xint = 0.101, RawFile="A11-12042"), aes(xintercept = xint), colour = "green", linewidth = 1) + 
  geom_vline(data = data.frame(xint = 0.0992, RawFile="A11-12043"), aes(xintercept = xint), colour = "green", linewidth = 1) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 16))

## Summary statistics for the Spectral angle per sample
Results %>%
  select(RawFile, SpectralAngle) %>%
  na.omit() %>%
  group_by(RawFile) %>%
  summarise(Avg_SA = mean(SpectralAngle, na.rm = TRUE),
            Median_SA = median(SpectralAngle, na.rm = TRUE),
            SD_SA = sd(SpectralAngle, na.rm = TRUE),
            Count_SA = n())

## Bootstrapping 95% Quantile intervals
### Seed
set.seed(1234)
### Sample A11-12042
Results %>%
  filter(RawFile == "A11-12042.mzML") %>%
  select(SpectralAngle) %>%
  na.omit() -> SA_12042
B=10000
t.boot<-c(1:B)
for(b in 1:B){
  x.boot<-sample(SA_12042$SpectralAngle, size=length(SA_12042$SpectralAngle),replace=T)
  t.boot[b]<-median(x.boot)
}
quantile(t.boot,probs=c(0.025,0.975))
### Sample A11-12043
Results %>%
  filter(RawFile == "A11-12043.mzML") %>%
  select(SpectralAngle) %>%
  na.omit() -> SA_12043
B=10000
t2.boot<-c(1:B)
for(b in 1:B){
  x2.boot<-sample(SA_12043$SpectralAngle, size=length(SA_12043$SpectralAngle),replace=T)
  t2.boot[b]<-median(x2.boot)
}
quantile(t2.boot,probs=c(0.025,0.975))

## Calculating the amount of detected peaks
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
### Summary table of amount of detected peaks
table(Results$NPeaks, Results$RawFile)

## Looking if the peaks are consecutive or not
Results["ConsecutivePeaks"] <- NA
for(index in 1:nrow(Results)){
  NPeaks = Results[[index, "NPeaks"]]
  if (NPeaks == 1) {
    Results[index, "ConsecutivePeaks"] <- "No"
  } else if(NPeaks == 2){
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]])){
      Results[index, "ConsecutivePeaks"] <- "Yes"
    } else {
      Results[index, "ConsecutivePeaks"] <- "No"
    }
  }  else if(NPeaks == 3){
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]])){
      Results[index, "ConsecutivePeaks"] <- "Yes"
    } else {
      Results[index, "ConsecutivePeaks"] <- "No"
    }
  } else if(NPeaks == 4){
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]])){
      Results[index, "ConsecutivePeaks"] <- "Yes"
    } else {
      Results[index, "ConsecutivePeaks"] <- "No"
    }
  } else if(NPeaks == 5){
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]]) & !is.na(Results[[index, "IsotopePeak5Intensity"]])){
      Results[index, "ConsecutivePeaks"] <- "Yes"
    } else {
      Results[index, "ConsecutivePeaks"] <- "No"
    }
  } else if(NPeaks == 6){
    if(!is.na(Results[[index, "IsotopePeak2Intensity"]]) & !is.na(Results[[index, "IsotopePeak3Intensity"]]) & !is.na(Results[[index, "IsotopePeak4Intensity"]]) & !is.na(Results[[index, "IsotopePeak5Intensity"]]) & !is.na(Results[[index, "IsotopePeak6Intensity"]])){
      Results[index, "ConsecutivePeaks"] <- "Yes"
    } else {
      Results[index, "ConsecutivePeaks"] <- "No"
    }
  }
}

### Table of amount of consecutive peaks
table(Results$NPeaks, Results$RawFile, Results$ConsecutivePeaks)

## Looking at the amount of unique peptide-modification-charge state combinations
data_unique <- unique(Results_final[ , c("PeptideSequence", "Modifications", "ChargeState")]) 

################################################################################
# Storing the final results into an Excel
Results %>%
  select(-NAs) -> Results_final
write_xlsx(Results_final, "Establishing A Comprehensive Workflow for Extracting MS1 Isotope Distributions in LC-MSMS Proteomics (Supporting Information 2).xlsx")

