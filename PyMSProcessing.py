"""
General script with functions for processing MS data in mzML format.
Written by Vilenne Frédérique
"""

# Importing functions
import os
import pywt
from pyopenms import *
import pandas as pd
from statsmodels import robust
from numpy import max, transpose
from numba import jit
import regex
from rpy2 import situation
os.environ["R_HOME"] = situation.get_r_home()
import rpy2.robjects as robjects
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import simplejson as json


class CustomJSONizer(json.JSONEncoder):
    def default(self, obj):
        return super().encode(bool(obj)) \
            if isinstance(obj, np.bool_) \
            else super().default(obj)

# Creation of functions
## Computing the mass
@jit
def mass_calculator(sequence: str, masses_type: str, cysteine_modified: int, methionine_modified: int):
    """
    Computes either the mono-isotopic mass or the average mass depending on the mass type input
    :param sequence: An amino-acid string of a peptide
    :param masses_type: The type of mass which has to be computed, either mono-isotopic or average.
    :param cysteine_modified: Amount of modified cysteines
    :param methionine_modified: Amount of modified methionines
    :return: The mass in Dalton, given a peptide string and mass type as a float.
    """
    # Setting dictionaries whom are required
    MonoIsotopicMass = {"A": 71.03711,
                        "R": 156.10111,
                        "N": 114.04293,
                        "D": 115.02694,
                        "C": 103.00919,
                        "E": 129.04259,
                        "Q": 128.05858,
                        "G": 57.02146,
                        "H": 137.05891,
                        "I": 113.08406,
                        "L": 113.08406,
                        "K": 128.09496,
                        "M": 131.04049,
                        "F": 147.06841,
                        "P": 97.05276,
                        "S": 87.03203,
                        "T": 101.04768,
                        "W": 186.07931,
                        "Y": 163.06333,
                        "V": 99.06841}

    AverageMass = {"A": 71.0788,
                   "R": 156.1875,
                   "N": 114.1038,
                   "D": 115.0886,
                   "C": 103.1388,
                   "E": 129.1155,
                   "Q": 128.1307,
                   "G": 57.0519,
                   "H": 137.1411,
                   "I": 113.1594,
                   "L": 113.1594,
                   "K": 128.1741,
                   "M": 131.1926,
                   "F": 147.1766,
                   "P": 97.1167,
                   "S": 87.0782,
                   "T": 101.1051,
                   "W": 186.2132,
                   "Y": 163.1760,
                   "V": 99.1326}

    Water = 18.0153
    Cysteine_mod = 57.0214  # Cysteine carbamidomethylation
    Methionine_mod = 15.9949    # Oxidation methionine

    sequence_mass = int(0)
    if masses_type == "MonoIsotopic":
        for AA in sequence:
            if AA in MonoIsotopicMass:
                sequence_mass += MonoIsotopicMass.get(AA)
            else:
                final_mass = sequence_mass + Water + Cysteine_mod * cysteine_modified + Methionine_mod * methionine_modified
                return final_mass
        final_mass = sequence_mass + Water + Cysteine_mod * cysteine_modified + Methionine_mod * methionine_modified
        return final_mass
    elif masses_type == "Average":
        for AA in sequence:
            if AA in AverageMass:
                sequence_mass += AverageMass.get(AA)
            else:
                final_mass = sequence_mass + Water + Cysteine_mod * cysteine_modified + Methionine_mod * methionine_modified
                return final_mass
        final_mass = sequence_mass + Water + Cysteine_mod * cysteine_modified + Methionine_mod * methionine_modified
        return final_mass


## Performing wavelet denoising
def wavelet_denoising(data, wavelet='db8'):
    """
    Perfoms Wavelet denoising according to https://doi.org/10.1002/pmic.200401261
    :param data: A dataframe consisting of the peak information extracted from a spectrum using getPeaks.
    :param wavelet: Dabauchies 8 wavelet
    :return: A smoothed spectrum
    """
    max_level = pywt.swt_max_level(len(data["Intensity"]))
    coefficients = pywt.swt(data=data["Intensity"], wavelet=wavelet, level=max_level, start_level=0, axis=-1, trim_approx=True, norm=True)
    mad = data[["Intensity"]].apply(robust.mad)
    mad = mad[0] * 6
    coefficients[1:] = (pywt.threshold(i, value=mad, mode='hard') for i in coefficients[1:])
    data["Intensity_denoised"] = pywt.iswt(coefficients, wavelet=wavelet)
    return data


## Performing wavelet denoising
def median_baseline_removal(data, noise_estimator=2):
    """
    Perfoms a median baseline removal based on the PepList algorithm from article https://doi.org/10.2174%2F138920209789177638
    :param data: A dataframe consisting of the peak information extracted from a spectrum using getPeaks.
    :param noise_estimator: How much the median is multiplied by for noise estimation
    :return: A baseline corrected spectrum
    """
    lower_limit = 360
    upper_limit = 407
    Jumps = 47
    while upper_limit <= 1300:
        data = data[data['MZ'].between(lower_limit, upper_limit)]
    return data


## Simple mass calculator
def mass_calculator_simple(sequence: str) -> float:
    """
    A more simplified version of the mass calculator which always computes the mono-isotopic mass.
    Only takes in a string as input so perfect for application on a Pandas dataframe.
    :param sequence: A string of amino-acids.
    :return: The monoisotopic mass as a float.
    """

    MonoIsotopicMass = {"A": 71.03711,
                        "R": 156.10111,
                        "N": 114.04293,
                        "D": 115.02694,
                        "C": 103.00919 + 57.02,  # Cysteine carbamidomethylation
                        "E": 129.04259,
                        "Q": 128.05858,
                        "G": 57.02146,
                        "H": 137.05891,
                        "I": 113.08406,
                        "L": 113.08406,
                        "K": 128.09496,
                        "M": 131.04049,
                        "F": 147.06841,
                        "P": 97.05276,
                        "S": 87.03203,
                        "T": 101.04768,
                        "W": 186.07931,
                        "Y": 163.06333,
                        "V": 99.06841}

    Water = 18.0153

    sequence_mass = int(0)
    for AA in sequence:
        if AA in MonoIsotopicMass:
            sequence_mass += MonoIsotopicMass.get(AA)
        else:
            final_mass = sequence_mass + Water
            return final_mass
    final_mass = sequence_mass + Water
    return final_mass


## Mass to Charge calculator
@jit
def mass_to_charge(mass: float, isotope: int, charge_state: int) -> float:
    """
    Computes the mass to charge value given the input.
    :param mass: The mass given, usually mono-isotopic mass in Dalton.
    :param isotope: Which isotope peak.
    :param charge_state: The charge state of the peptide.
    :return: The M/Z value of the peptide.
    """
    Hydrogen = 1.009
    mz = (mass + isotope + Hydrogen * charge_state) / charge_state
    return mz


@jit
def mass_to_charge_observed_mz(mz_observed: float, isotope: int, charge_state: int) -> float:
    """
    Computes the mass to charge value given the input.
    :param mz_observed: The observed mz
    :param isotope: Which isotope peak.
    :param charge_state: The charge state of the peptide.
    :return: The M/Z value of the peptide.
    """
    Hydrogen = 1.009
    mz = mz_observed + isotope / charge_state
    return mz


## Extracted Ion Chromatogram
def xic(data, mz: float, tolerance: int):
    """
    Generates an Extracted Ion Chromatogram
    :param data: Data in mzML format loaded in using PyOpenMS
    :param mz: A value for MZ
    :param tolerance: A tolerance for ppm, usually 1 or 5
    :return: Returns a Pandas dataframe with the Retention times and maximum intensities within the MZ tolerance
    """
    xic_df = pd.DataFrame()
    spectral_id = []
    retention_times = []
    intensities = []
    tolerance_window = mz / (1000000 / tolerance)
    lower_limit = mz - tolerance_window
    upper_limit = mz + tolerance_window
    for spec in data:
        spectral_id.append(spec.getNativeID())
        retention_times.append(spec.getRT())
        array_temp = spec.get_peaks()
        array_temp = transpose(array_temp)
        array_temp2 = array_temp[np.where((array_temp[:, 0] >= lower_limit) & (array_temp[:, 0] <= upper_limit))]
        try:
            max_intensity = max(array_temp2[:, 1])
        except ValueError:
            max_intensity = 0
            intensities.append(max_intensity)
        else:
            intensities.append(max_intensity)
    xic_df["SpectraID"] = spectral_id
    xic_df["Retention Time"] = retention_times
    xic_df["Intensity"] = intensities
    return xic_df


## Extract isotopes
@jit
def isotope_extractor(peaks, storage_list, mz, charge, tolerance):
    """
    Generates a Total Ion Chromatogram
    :param peaks: Peak list with intensities
    :param storage_list: List to add isotope to
    :param mz: MZ Value to look at
    :param charge: Charge state of peptide
    :param tolerance: tolerance in ppm
    :return: Returns the list with additional isotope distribution
    """
    IsotopePeaks = [0, 1, 2, 3, 4, 5]
    for isotope in IsotopePeaks:
        IsotopicPeakMZ = mass_to_charge_observed_mz(mz, isotope=isotope, charge_state=charge)

        ## 5 ppm tolerance
        tolerance_window = IsotopicPeakMZ / (1000000 / tolerance)
        lower_limit, upper_limit = IsotopicPeakMZ - tolerance_window, IsotopicPeakMZ + tolerance_window
        array_temp2 = peaks[np.where((peaks[:, 0] >= lower_limit) & (peaks[:, 0] <= upper_limit))]

        if len(array_temp2) == 0:
            MZValue, IntensityValue = np.nan, np.nan
            storage_list.extend([MZValue, IntensityValue])

        else:
            Peak_information = array_temp2[np.where(array_temp2[:, 1] == max(array_temp2[:, 1]))]
            MZValue, IntensityValue = Peak_information[0, 0], Peak_information[0, 1]
            storage_list.extend([MZValue, IntensityValue])

        # try:
        #     Peak_information = array_temp2[np.where(array_temp2[:, 1] == max(array_temp2[:, 1]))]
        # except ValueError:
        #     MZValue, IntensityValue = np.nan, np.nan
        #     storage_list.extend([MZValue, IntensityValue])
        # else:
        #     MZValue, IntensityValue = Peak_information[0, 0], Peak_information[0, 1]
        #     storage_list.extend([MZValue, IntensityValue])
    return storage_list


## Total Ion Chromatogram
def tic(data):
    """
    Generates a Total Ion Chromatogram
    :param data: Data in mzML format loaded in using PyOpenMS
    :return: Returns a Pandas dataframe with the Retention times and total ion current at the retention time
    """
    tic_df = pd.DataFrame()
    retention_times = []
    tic_value = []
    for spectra in data:
        retention_times.append(spectra.getRT())
        tic_value.append(spectra.calculateTIC())
    tic_df["Retention Time"] = retention_times
    tic_df["TIC"] = tic_value
    return tic_df


## Baseline Ion Chromatogram
def bic(data):
    """
    Generates the data for a baseline ion chromatogram
    :param data: Data in mzML format loaded in using PyOpenMS
    :return: Returns a Pandas dataframe with the needed data
    """
    bic_df = pd.DataFrame()
    retention_times = []
    intensity = []
    mz = []
    for spectra in data:
        retention_times.append(spectra.getRT())
        array_temp = spectra.get_peaks()
        array_temp = np.transpose(array_temp)
        df_temp = pd.DataFrame(data=array_temp, columns=['MZ', 'Intensity'])
        # Get the highest intensity
        df_temp = df_temp[df_temp.Intensity == df_temp.Intensity.max()]
        intensity.append(df_temp["Intensity"].values[0])
        mz.append(df_temp["MZ"].values[0])
    bic_df["Retention Time"] = retention_times
    bic_df["Intensity"] = intensity
    bic_df["MZ"] = mz
    return bic_df


## Heatmap construction
def total_heatmap(data):
    """
    Generates the data for a heatmap of all observations
    :param data: Data in mzML format loaded in using PyOpenMS
    :return: Returns a Pandas dataframe with the needed data, looks shite though
    """
    heatmap_df = pd.DataFrame()
    retention_times = []
    total_ion_counts = []
    intensity = []
    mz = []
    for spectra in data:
        retention_time = spectra.getRT()
        total_ion_count = spectra.calculateTIC()
        array_temp = spectra.get_peaks()
        array_temp = np.transpose(array_temp)
        df_temp = pd.DataFrame(data=array_temp, columns=['MZ', 'Intensity'])
        length = len(df_temp.index)
        retention_times.extend([retention_time for i in range(length)])
        total_ion_counts.extend([total_ion_count for i in range(length)])
        intensity.extend(df_temp["Intensity"])
        mz.extend(df_temp["MZ"])
    heatmap_df["Retention Time"] = retention_times
    heatmap_df["TIC"] = total_ion_counts
    heatmap_df["Intensity"] = intensity
    heatmap_df["MZ"] = mz
    return heatmap_df


## MS level filter
def ms_level_filter(input_directory: str, ms_level: int):
    """
    Filters the mzML file for the requested ms-level
    :param input_directory: A string to the location of the mzML file
    :param ms_level: Integer of the ms-level, usually 1 or 2
    :return: Filtered mzML data
    """
    data = MSExperiment()
    MzMLFile().load(input_directory, data)
    spec = []
    for s in data.getSpectra():
        if s.getMSLevel() == ms_level:
            spec.append(s)
    data.setSpectra(spec)
    MzMLFile().store(input_directory, data)
    return


# Trypsin cleaver
def trypsin_digest(sequence):
    Start_index = 0
    Fragments = []
    for aa in range(0, len(sequence)):
        if aa == len(sequence) - 1:
            Fragments += [sequence[Start_index:aa + 1]]
        elif sequence[aa] == "R" and sequence[aa + 1] != "P":
            Fragments += [sequence[Start_index:aa + 1]]
            Start_index = aa + 1
        elif sequence[aa] == "K" and sequence[aa + 1] != "P":
            Fragments += [sequence[Start_index:aa + 1]]
            Start_index = aa + 1
    return Fragments


def UPS_concentration_lookup(proteinID: str):
    Concentrations = {
        "P02768ups|ALBU_HUMAN_UPS Serum albumin (Chain 26-609) - Homo sapiens (Human)": 50,
        "Q15843ups|NEDD8_HUMAN_UPS NEDD8 (Chain 1-81) - Homo sapiens (Human)": 0.5,
        "P01112ups|RASH_HUMAN_UPS GTPase HRas (Chain 1-189) - Homo sapiens (Human)": 0.005,
        "P04040ups|CATA_HUMAN_UPS Catalase (Chain 2-527) - Homo sapiens (Human)": 5,
        "P02753ups|RETBP_HUMAN_UPS Retinol-binding protein 4 (Chain 19-201) - Homo sapiens (Human)": 0.5,
        "P06396ups|GELS_HUMAN_UPS Gelsolin (Chain 28-782) - Homo sapiens (Human)": 0.005,
        "P01375ups|TNFA_HUMAN_UPS Tumor necrosis factor, soluble form (Chain 77-233) - Homo sapiens (Human)": 0.0005,
        "P15559ups|NQO1_HUMAN_UPS NAD(P)H dehydrogenase [quinone] 1 (Chain 2-274) - Homo sapiens (Human)": 5,
        "Q06830ups|PRDX1_HUMAN_UPS Peroxiredoxin 1 (Chain 2-199) - Homo sapiens (Human)": 5,
        "P00167ups|CYB5_HUMAN_UPS Cytochrome b5 (Chain 1-134, N-terminal His tag) - Homo sapiens (Human)": 5,
        "P06732ups|KCRM_HUMAN_UPS Creatine kinase M-type (Chain 1-381) - Homo sapiens (Human)": 0.5,
        "P02741ups|CRP_HUMAN_UPS C-reactive protein (Chain 19-224) - Homo sapiens (Human)": 0.0005,
        "P61626ups|LYSC_HUMAN_UPS Lysozyme C (Chain 19-148) - Homo sapiens (Human)": 0.5,
        "P16083ups|NQO2_HUMAN_UPS Ribosyldihydronicotinamide dehydrogenase [quinone] (Chain 2-231) - Homo sapiens (Human)": 0.5,
        "P10145ups|IL8_HUMAN_UPS Interleukin-8, IL-8 (Chain 28-99) - Homo sapiens (Human)": 0.0005,
        "P61769ups|B2MG_HUMAN_UPS Beta-2-microglobulin (Chain 21-119) - Homo sapiens (Human)": 0.05,
        "P02144ups|MYG_HUMAN_UPS Myoglobin (Chain 2-154) - Homo sapiens (Human)": 5,
        "P08263ups|GSTA1_HUMAN_UPS Glutathione S-transferase A1 (Chain 2-222) - Homo sapiens (Human)": 0.05,
        "P55957ups|BID_HUMAN_UPS BH3-interacting domain death agonist (Chain 1-195) - Homo sapiens (Human)": 0.05,
        "P00709ups|LALBA_HUMAN_UPS Alpha-lactalbumin (Chain 20-142) - Homo sapiens (Human)": 0.05,
        "P69905ups|HBA_HUMAN_UPS Hemoglobin subunit alpha (Chain 2-142) - Homo sapiens (Human)": 50,
        "P01344ups|IGF2_HUMAN_UPS Insulin-like growth factor II (Chain 25-91) - Homo sapiens (Human)": 0.05,
        "P00915ups|CAH1_HUMAN_UPS Carbonic anhydrase 1 (Chain 2-261) - Homo sapiens (Human)": 50,
        "P08758ups|ANXA5_HUMAN_UPS Annexin A5 (Chain 2-320) - Homo sapiens (Human)": 0.0005,
        "P00918ups|CAH2_HUMAN_UPS Carbonic anhydrase 2 (Chain 2-260) - Homo sapiens (Human)": 50,
        "P62937ups|PPIA_HUMAN_UPS Peptidyl-prolyl cis-trans isomerase A (Chain 1-165, N terminal His tag)- Homo sapiens (Human)": 5,
        "P00441ups|SODC_HUMAN_UPS Superoxide dismutase [Cu-Zn] (Chain 2-154) - Homo sapiens (Human)": 0.0005,
        "P05413ups|FABPH_HUMAN_UPS Fatty acid-binding protein, heart (Chain 2-133) - Homo sapiens (Human)": 0.0005,
        "O00762ups|UBE2C_HUMAN_UPS Ubiquitin-conjugating enzyme E2 C (Chain 1-179, N-terminal His tag)- Homo sapiens (Human)": 0.005,
        "P41159ups|LEP_HUMAN_UPS Leptin (Chain 22-167) - Homo sapiens (Human)": 50,
        "P02788ups|TRFL_HUMAN_UPS Lactotransferrin (Chain 20-710) - Homo sapiens (Human)": 0.0005,
        "P09211ups|GSTP1_HUMAN_UPS Glutathione S-transferase P (Chain 2-210) - Homo sapiens (Human)": 0.05,
        "P01031ups|CO5_HUMAN_UPS Complement C5 (C5a anaphylatoxin) (Chain 678-751) - Homo sapiens (Human)": 50,
        "P10636-8ups|TAU_HUMAN_UPS Microtubule-associated protein tau {Isoform Tau-F (Tau-4)} (Chain 2-441) - Homo sapiens (Human)": 0.0005,
        "P01133ups|EGF_HUMAN_UPS Pro-Epidermal growth factor (EGF) (Chain 971-1023) - Homo sapiens (Human)": 5,
        "P63165ups|SUMO1_HUMAN_UPS Small ubiquitin-related modifier 1 (Chain 1-97, N-terminal GST tag) - Homo sapiens (Human)": 5,
        "P51965ups|UB2E1_HUMAN_UPS Ubiquitin-conjugating enzyme E2 E1 (Chain 1-193, N terminal His tag)- Homo sapiens (Human)": 0.005,
        "P12081ups|SYHC_HUMAN_UPS Histidyl-tRNA synthetase, cytoplasmic (Chain 1-509, C terminal His tag) - Homo sapiens (Human)": 0.5,
        "P99999ups|CYC_HUMAN_UPS Cytochrome c (Chain 2-105) - Homo sapiens (Human)": 0.005,
        "P02787ups|TRFE_HUMAN_UPS Serotransferrin (Chain 20-698) - Homo sapiens (Human)": 0.005,
        "P10599ups|THIO_HUMAN_UPS Thioredoxin (Chain 2-105, N-terminal His tag)- Homo sapiens (Human)": 0.05,
        "P62988ups|UBIQ_HUMAN_UPS Ubiquitin (Chain 1-76, N-terminal His tag) - Homo sapiens (Human)": 50,
        "P01008ups|ANT3_HUMAN_UPS Antithrombin-III (Chain 33-464) - Homo sapiens (Human)": 0.05,
        "P01127ups|PDGFB_HUMAN_UPS Platelet-derived growth factor B chain (Chain 82-190) - Homo sapiens (Human)": 0.05,
        "P68871ups|HBB_HUMAN_UPS Hemoglobin subunit beta (Chain 2-147) - Homo sapiens (Human)": 50,
        "P63279ups|UBC9_HUMAN_UPS SUMO-conjugating enzyme UBC9 (Chain 1-158) - Homo sapiens (Human)": 0.5,
        "O76070ups|SYUG_HUMAN_UPS Gamma-synuclein (Chain 1-127) - Homo sapiens (Human)": 0.05,
        "P01579ups|IFNG_HUMAN_UPS Interferon Gamma (Chain 23-166) - Homo sapiens (Human)": 0.005
    }
    Concentration = Concentrations[proteinID]
    return Concentration


def excel_MS1_Isotope_Distributions(raw_file, input_list, output_file):
    # Storage
    Storage = []

    # Opening the raw MS data
    Data = MSExperiment()
    MzMLFile().load(raw_file, Data)

    # Reading the input list
    input_data = pd.read_excel(input_list)

    # Check if columns are present in the input data
    if 'Peptide_Sequence' in input_data.columns:
        print("Peptide_Sequence column is present")
    else:
        print("Peptide_Sequence column is absent, algorithm will abort")
        return
    if 'Charge_State' in input_data.columns:
        print("Charge_State column is present")
    else:
        print("Charge_State column is absent, algorithm will abort")
        return
    if 'Modifications' in input_data.columns:
        print("Modifications column is present")
    else:
        print("Modifications column is absent, algorithm will abort")
        return
    if 'RetentionTime' in input_data.columns:
        print("RetentionTime column is present")
    else:
        print("RetentionTime column is absent, algorithm will abort")
        return
    if 'RetentionTimeWindowBefore' in input_data.columns:
        print("RetentionTimeWindowBefore column is present")
    else:
        print("RetentionTimeWindowBefore column is absent, algorithm will abort")
        return
    if 'RetentionTimeWindowAfter' in input_data.columns:
        print("RetentionTimeWindowAfter column is present")
    else:
        print("RetentionTimeWindowAfter column is absent, algorithm will abort")
        return


    for index in range(len(input_data)):
        # Extract metadata
        Peptide = input_data["Peptide_Sequence"].iloc[index]
        Charge = input_data["Charge_State"].iloc[index]
        Modifications = input_data["Modifications"].iloc[index]
        RetentionTime = input_data["RetentionTime"].iloc[index]
        RetentionTimeWindowBefore = input_data["RetentionTimeWindowBefore"].iloc[index]
        RetentionTimeWindowAfter = input_data["RetentionTimeWindowAfter"].iloc[index]

        PeptideLength = len(Peptide)

        # Amount of Modifications
        if pd.isnull(Modifications):
            CysteineMod, MethionineMod = 0, 0
            TheoreticalMonoIsotopicMass, TheoreticalAverageMass = mass_calculator(sequence=Peptide,
                                                                                  masses_type="MonoIsotopic",
                                                                                  cysteine_modified=CysteineMod,
                                                                                  methionine_modified=MethionineMod), mass_calculator(
                sequence=Peptide, masses_type="Average", cysteine_modified=CysteineMod,
                methionine_modified=MethionineMod)

        else:
            CysteineMod, MethionineMod = len(regex.findall(pattern=r"C\(57.0214\)", string=Modifications)), len(
                regex.findall(pattern=r"M\(15.9949\)", string=Modifications))
            TheoreticalMonoIsotopicMass, TheoreticalAverageMass = mass_calculator(sequence=Peptide,
                                                                                  masses_type="MonoIsotopic",
                                                                                  cysteine_modified=CysteineMod,
                                                                                  methionine_modified=MethionineMod), mass_calculator(
                sequence=Peptide, masses_type="Average", cysteine_modified=CysteineMod,
                methionine_modified=MethionineMod)

        TheoreticalMZ = mass_to_charge(mass=TheoreticalMonoIsotopicMass, isotope=0, charge_state=Charge)

        if 'ErrorTolerance' in input_data.columns:
            print("ErrorTolerance column is present")
            ErrorTolerance = input_data["ErrorTolerance"].iloc[index]
        else:
            print("ErrorTolerance column is absent, a default value of 5ppm will be used")
            ErrorTolerance = 5
        if 'MZ' in input_data.columns:
            print("MZ column is present")
            MZ = input_data["MZ"].iloc[index]
        else:
            print("MZ column is absent, the theoretical MZ values will be used")
            MZ = TheoreticalMZ

        # Extracting Ion Chromatogram
        XIC_results = xic(data=Data, mz=MZ, tolerance=ErrorTolerance)

        # Extracting data from XIC
        ## Filtering XIC for only relevant spectra
        XIC_results = XIC_results[XIC_results['Intensity'] > 0]
        XIC_results = XIC_results[XIC_results['Retention Time'].between(RetentionTime - RetentionTimeWindowBefore, RetentionTime + RetentionTimeWindowAfter)]

        ## Amount of spectra
        NSpectra = len(XIC_results.index)
        if NSpectra > 0:
            # Continuing for isotope distributions

            ## Extracting all spectra ID
            SpectraID = XIC_results["SpectraID"].tolist()

            # Iterate over all spectra from XIC
            for SpectrumID in SpectraID:
                ## Prepping list for storage
                Storage_temp = []

                ## Looping over all spectra in RAW-file
                for Spectrum in Data:
                    Spectrum_ID_Temp = Spectrum.getNativeID()
                    if Spectrum_ID_Temp == SpectrumID:
                        ## Get final summary stats
                        RT, TIC = Spectrum.getRT(), Spectrum.calculateTIC()

                        ## Appending to list before iterating over peaks
                        Storage_temp.extend(
                            [raw_file, SpectrumID, Peptide, Modifications, PeptideLength, TheoreticalMonoIsotopicMass,
                             TheoreticalAverageMass, Charge, TheoreticalMZ, MZ, RT, TIC])

                        # Acquiring peak information
                        ## Acquiring the peak information from spectrum
                        array_temp = np.transpose(Spectrum.get_peaks())

                        # Get isotopes
                        Storage_temp = isotope_extractor(peaks=array_temp, storage_list=Storage_temp, mz=MZ, charge=Charge, tolerance=ErrorTolerance)

                        # Add to final storage list
                        Storage.append(Storage_temp)

    # Generate the Pandas dataframe
    Results = pd.DataFrame(Storage, columns=['RawFile',
                                             'Spectrum',
                                             'PeptideSequence',
                                             'Modifications',
                                             'PeptideLength',
                                             'TheoreticalMonoIsotopicMass',
                                             'TheoreticalAverageMass',
                                             'ChargeState',
                                             'TheoreticalMZ',
                                             'ObservedMZ',
                                             'RetentionTime',
                                             'TIC',
                                             'IsotopePeak1MZ',
                                             'IsotopePeak1Intensity',
                                             'IsotopePeak2MZ',
                                             'IsotopePeak2Intensity',
                                             'IsotopePeak3MZ',
                                             'IsotopePeak3Intensity',
                                             'IsotopePeak4MZ',
                                             'IsotopePeak4Intensity',
                                             'IsotopePeak5MZ',
                                             'IsotopePeak5Intensity',
                                             'IsotopePeak6MZ',
                                             'IsotopePeak6Intensity'])

    # Setting Modifications correct
    Results["Modifications"] = Results["Modifications"].fillna("None")

    # Start processing in R for BRAIN
    print("Further processing in R, using BRAIN")
    # Defining the R script and loading the instance in Python
    r = robjects.r
    r['source']('MS1ProcessingR.R')
    # Loading the function
    MS1_Isotope_Distribution_Function_R = robjects.globalenv['MS1_Isotope_Distribution_Function']
    # converting it into r object for passing into r function
    with (ro.default_converter + pandas2ri.converter).context():
        Results_R = ro.conversion.get_conversion().py2rpy(Results)
    # Invoking the R function and getting the result
    Results_R = MS1_Isotope_Distribution_Function_R(Results_R)
    # Converting it back to a pandas dataframe.
    with (ro.default_converter + pandas2ri.converter).context():
        Results_Py= ro.conversion.get_conversion().rpy2py(Results_R)

    # Write in final results
    Results_Py.to_excel(f"{output_file}.xlsx", index=False)

    # Done processing
    print("Done Processing")

    return


# JSON Function
def JSON_MS1_Isotope_Distributions(raw_file, input_list, output_file):
    # Start out by creating the Excel file
    excel_MS1_Isotope_Distributions(raw_file=raw_file, input_list=input_list, output_file=output_file)

    # Read the Excel
    Data = pd.read_excel(f"{output_file}.xlsx")

    # Setting up a dictionary to store the results in
    JSON_Dict = {}

    # Acquire all unique peptides
    Peptides = list(set(Data["PeptideSequence"].tolist()))

    # Iterate over all Peptides
    for Peptide in Peptides:
        # Filter the dataframe for the corresponding peptide
        Data_Peptide = Data[Data["PeptideSequence"] == Peptide]

        # Extract all charge States
        Charge_States = list(set(Data_Peptide["ChargeState"].tolist()))

        # Extract parameters from data frame
        PeptideLength = Data_Peptide["PeptideLength"].iloc[0].item()

        # Add information to Dictionary
        JSON_Dict[Peptide] = {"Peptide Metadata": {"Peptide Length": PeptideLength}}

        # Temporary dictionary to store results in
        Temp_JSON_Dict_CS = {}

        # Iterate over all Charge States
        for Charge_State in Charge_States:
            # Temporary dictionary to store results in
            Temp_JSON_Dict_Modifications = {}

            # Charge State string for dictionary
            Charge_State_Dict = f"Charge State +{Charge_State}"

            # Filter the dataframe for the corresponding Charge State
            Data_Peptide_CS = Data_Peptide[Data_Peptide["ChargeState"] == Charge_State]

            # Extract all modifications
            Modifications = list(set(Data_Peptide_CS["Modifications"].tolist()))

            # Iterate over all modifications
            for Modification in Modifications:
                # Filter the dataframe for the corresponding Modification, having to deal with NA (no modification) and modifications
                Temp_JSON_Dict_Modifications_To_append = {}

                if pd.isna(Modification):
                    # Filter the dataframe for the corresponding Modification in case of NA (Modification Absent)
                    Data_Temp_CS_Mod = Data_Peptide_CS[Data_Peptide_CS["Modifications"].isna()]

                    # Modification string for dictionary
                    Modification_Dict = f"Modification: None"

                    # Extract information
                    TheoreticalMonoIsotopicMass = Data_Temp_CS_Mod["TheoreticalMonoIsotopicMass"].iloc[0].item()
                    TheoreticalAverageMass = Data_Temp_CS_Mod["TheoreticalAverageMass"].iloc[0].item()
                    TheoreticalMZ = Data_Temp_CS_Mod["TheoreticalMZ"].iloc[0].item()
                    ObservedMZ = Data_Temp_CS_Mod["ObservedMZ"].iloc[0].item()
                    Carbons = Data_Temp_CS_Mod["Carbons"].iloc[0].item()
                    Hydrogens = Data_Temp_CS_Mod["Hydrogens"].iloc[0].item()
                    Oxygens = Data_Temp_CS_Mod["Oxygens"].iloc[0].item()
                    Nitrogens = Data_Temp_CS_Mod["Nitrogens"].iloc[0].item()
                    Sulphurs = Data_Temp_CS_Mod["Sulphurs"].iloc[0].item()

                    # Temporary dictionary to store results in
                    Temp_JSON_Dict_Modifications_To_append = {Modification_Dict: {
                        "Ion Metadata": {
                            "TheoreticalMonoIsotopicMass": TheoreticalMonoIsotopicMass,
                            "TheoreticalAverageMass": TheoreticalAverageMass,
                            "TheoreticalMZ": TheoreticalMZ,
                            "ObservedMZ": ObservedMZ,
                            "Carbons": Carbons,
                            "Hydrogens": Hydrogens,
                            "Oxygens": Oxygens,
                            "Nitrogens": Nitrogens,
                            "Sulphurs": Sulphurs,
                        }
                    }
                    }

                    # Start iterating over dataframe to acquire information
                    for index in range(len(Data_Temp_CS_Mod)):
                        # Extract all information
                        Spectrum = Data_Temp_CS_Mod["Spectrum"].iloc[index]
                        RetentionTime = Data_Temp_CS_Mod["RetentionTime"].iloc[index].item()
                        TIC = Data_Temp_CS_Mod["TIC"].iloc[index].item()
                        Distribution = Data_Temp_CS_Mod["Distribution"].iloc[index]
                        SpectralAngle = Data_Temp_CS_Mod["SpectralAngle"].iloc[index].item()
                        NPeaks = Data_Temp_CS_Mod["NPeaks"].iloc[index].item()
                        Consecutive_Peaks = Data_Temp_CS_Mod["ConsecutivePeaks"].iloc[index]

                        # Store isotope distribution
                        Isotope = {
                            "MZ": {"IsotopePeak1MZ": Data_Temp_CS_Mod["IsotopePeak1MZ"].iloc[index].item(),
                                   "IsotopePeak2MZ": Data_Temp_CS_Mod["IsotopePeak2MZ"].iloc[index].item(),
                                   "IsotopePeak3MZ": Data_Temp_CS_Mod["IsotopePeak3MZ"].iloc[index].item(),
                                   "IsotopePeak4MZ": Data_Temp_CS_Mod["IsotopePeak4MZ"].iloc[index].item(),
                                   "IsotopePeak5MZ": Data_Temp_CS_Mod["IsotopePeak5MZ"].iloc[index].item(),
                                   "IsotopePeak6MZ": Data_Temp_CS_Mod["IsotopePeak6MZ"].iloc[index].item()},
                            "Intensities": {
                                "IsotopePeak1Intensity": Data_Temp_CS_Mod["IsotopePeak1Intensity"].iloc[index].item(),
                                "IsotopePeak2Intensity": Data_Temp_CS_Mod["IsotopePeak2Intensity"].iloc[index].item(),
                                "IsotopePeak3Intensity": Data_Temp_CS_Mod["IsotopePeak3Intensity"].iloc[index].item(),
                                "IsotopePeak4Intensity": Data_Temp_CS_Mod["IsotopePeak4Intensity"].iloc[index].item(),
                                "IsotopePeak5Intensity": Data_Temp_CS_Mod["IsotopePeak5Intensity"].iloc[index].item(),
                                "IsotopePeak6Intensity": Data_Temp_CS_Mod["IsotopePeak6Intensity"].iloc[index].item()}}

                        # Store BRAIN Results
                        BRAIN = {"BRAINRelativeIsotopePeak1Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak1Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak2Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak2Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak3Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak3Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak4Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak4Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak5Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak5Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak6Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak6Intensity"].iloc[index].item()}

                        Temp_JSON_Dict_Modifications_To_append[Modification_Dict][Spectrum] = {
                            "RetentionTime": RetentionTime,
                            "TIC": TIC,
                            "Distribution": Distribution,
                            "NPeaks": NPeaks,
                            "ConsecutivePeaks": Consecutive_Peaks,
                            "SpectralAngle": SpectralAngle,
                            "IsotopeDistribution": Isotope,
                            "BRAINDistribution": BRAIN
                        }


                else:
                    # Filter the dataframe for the corresponding Modification in case of non-NA (Modification present)
                    Data_Temp_CS_Mod = Data_Peptide_CS[Data_Peptide_CS["Modifications"] == Modification]

                    # Modification string for dictionary
                    Modification_Dict = f"Modification: {Modification}"

                    # Extract information
                    TheoreticalMonoIsotopicMass = Data_Temp_CS_Mod["TheoreticalMonoIsotopicMass"].iloc[0].item()
                    TheoreticalAverageMass = Data_Temp_CS_Mod["TheoreticalAverageMass"].iloc[0].item()
                    TheoreticalMZ = Data_Temp_CS_Mod["TheoreticalMZ"].iloc[0].item()
                    ObservedMZ = Data_Temp_CS_Mod["ObservedMZ"].iloc[0].item()
                    Carbons = Data_Temp_CS_Mod["Carbons"].iloc[0].item()
                    Hydrogens = Data_Temp_CS_Mod["Hydrogens"].iloc[0].item()
                    Oxygens = Data_Temp_CS_Mod["Oxygens"].iloc[0].item()
                    Nitrogens = Data_Temp_CS_Mod["Nitrogens"].iloc[0].item()
                    Sulphurs = Data_Temp_CS_Mod["Sulphurs"].iloc[0].item()

                    # Temporary dictionary to store results in
                    Temp_JSON_Dict_Modifications_To_append = {Modification_Dict: {
                        "Ion Metadata": {
                            "TheoreticalMonoIsotopicMass": TheoreticalMonoIsotopicMass,
                            "TheoreticalAverageMass": TheoreticalAverageMass,
                            "TheoreticalMZ": TheoreticalMZ,
                            "ObservedMZ": ObservedMZ,
                            "Carbons": Carbons,
                            "Hydrogens": Hydrogens,
                            "Oxygens": Oxygens,
                            "Nitrogens": Nitrogens,
                            "Sulphurs": Sulphurs,
                        }
                    }
                    }

                    # Start iterating over dataframe to acquire information
                    for index in range(len(Data_Temp_CS_Mod)):
                        # Extract all information
                        Spectrum = Data_Temp_CS_Mod["Spectrum"].iloc[index]
                        RetentionTime = Data_Temp_CS_Mod["RetentionTime"].iloc[index].item()
                        TIC = Data_Temp_CS_Mod["TIC"].iloc[index].item()
                        Distribution = Data_Temp_CS_Mod["Distribution"].iloc[index]
                        SpectralAngle = Data_Temp_CS_Mod["SpectralAngle"].iloc[index].item()
                        NPeaks = Data_Temp_CS_Mod["NPeaks"].iloc[index].item()
                        Consecutive_Peaks = Data_Temp_CS_Mod["ConsecutivePeaks"].iloc[index]

                        # Store isotope distribution
                        Isotope = {
                            "MZ": {"IsotopePeak1MZ": Data_Temp_CS_Mod["IsotopePeak1MZ"].iloc[index].item(),
                                   "IsotopePeak2MZ": Data_Temp_CS_Mod["IsotopePeak2MZ"].iloc[index].item(),
                                   "IsotopePeak3MZ": Data_Temp_CS_Mod["IsotopePeak3MZ"].iloc[index].item(),
                                   "IsotopePeak4MZ": Data_Temp_CS_Mod["IsotopePeak4MZ"].iloc[index].item(),
                                   "IsotopePeak5MZ": Data_Temp_CS_Mod["IsotopePeak5MZ"].iloc[index].item(),
                                   "IsotopePeak6MZ": Data_Temp_CS_Mod["IsotopePeak6MZ"].iloc[index].item()},
                            "Intensities": {
                                "IsotopePeak1Intensity": Data_Temp_CS_Mod["IsotopePeak1Intensity"].iloc[index].item(),
                                "IsotopePeak2Intensity": Data_Temp_CS_Mod["IsotopePeak2Intensity"].iloc[index].item(),
                                "IsotopePeak3Intensity": Data_Temp_CS_Mod["IsotopePeak3Intensity"].iloc[index].item(),
                                "IsotopePeak4Intensity": Data_Temp_CS_Mod["IsotopePeak4Intensity"].iloc[index].item(),
                                "IsotopePeak5Intensity": Data_Temp_CS_Mod["IsotopePeak5Intensity"].iloc[index].item(),
                                "IsotopePeak6Intensity": Data_Temp_CS_Mod["IsotopePeak6Intensity"].iloc[index].item()}}

                        # Store BRAIN Results
                        BRAIN = {"BRAINRelativeIsotopePeak1Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak1Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak2Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak2Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak3Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak3Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak4Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak4Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak5Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak5Intensity"].iloc[index].item(),
                                 "BRAINRelativeIsotopePeak6Intensity":
                                     Data_Temp_CS_Mod["BRAINRelativeIsotopePeak6Intensity"].iloc[index].item()}

                        Temp_JSON_Dict_Modifications_To_append[Modification_Dict][Spectrum] = {
                            "RetentionTime": RetentionTime,
                            "TIC": TIC,
                            "Distribution": Distribution,
                            "NPeaks": NPeaks,
                            "ConsecutivePeaks": Consecutive_Peaks,
                            "SpectralAngle": SpectralAngle,
                            "IsotopeDistribution": Isotope,
                            "BRAINDistribution": BRAIN
                        }

                Temp_JSON_Dict_Modifications.update(Temp_JSON_Dict_Modifications_To_append)

            Temp_JSON_Dict_CS[Charge_State_Dict] = Temp_JSON_Dict_Modifications

        # Have to add to the final JSON after each iteration
        JSON_Dict[Peptide].update(Temp_JSON_Dict_CS)

    # Serializing json
    Final_Dict = {"Peptides": JSON_Dict}
    json_object = json.dumps(Final_Dict, indent=6, ignore_nan=True, cls=CustomJSONizer)

    # Writing to sample.json
    with open(f"{output_file}.json", "w") as outfile:
        outfile.write(json_object)
