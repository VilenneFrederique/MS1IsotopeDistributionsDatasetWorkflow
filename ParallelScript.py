# Modules
import pandas as pd
import numpy as np
import regex
from pyopenms import *
from tqdm import tqdm
import os
import PyMSProcessing
from Bio import SeqIO


## Iterate over all files
def processor(file):
    ### Extracting parameter raw file
    raw_file = file

    # Storage
    Storage = []

    ### Opening the data
    Data = MSExperiment()
    MzMLFile().load(f"G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MS1 Spectra/{raw_file}", Data)

    ### Results MSFragger
    Sample = regex.findall(pattern=r"A11-\d*", string=file)[0]
    MSFraggerResults = pd.read_table(f"G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/MSFragger/{Sample}/psm.tsv")

    ### Create numpy array
    MSFraggerResults_list = MSFraggerResults.values.tolist()

    for rows in tqdm(MSFraggerResults_list):

        # Acquire several statistics
        Peptide, PeptideLength, Charge, ObservedMZ, PeptideProphetScore, Modifications, ProteinID = rows[2], rows[6], rows[7], rows[11], rows[21], rows[27], rows[31]

        # RT
        RT_Peptide = rows[8]

        # Get protein concentration
        Concentration = PyMSProcessing.UPS_concentration_lookup(proteinID=ProteinID)

        # Amount of Modifications
        if pd.isnull(Modifications):
            CysteineMod, MethionineMod = 0, 0
            TheoreticalMonoIsotopicMass, TheoreticalAverageMass = PyMSProcessing.mass_calculator(sequence=Peptide,
                                                                                                 masses_type="MonoIsotopic",
                                                                                                 cysteine_modified=CysteineMod,
                                                                                                 methionine_modified=MethionineMod), PyMSProcessing.mass_calculator(
                sequence=Peptide, masses_type="Average", cysteine_modified=CysteineMod,
                methionine_modified=MethionineMod)

        else:
            CysteineMod, MethionineMod = len(regex.findall(pattern=r"C\(57.0214\)", string=Modifications)), len(
                regex.findall(pattern=r"M\(15.9949\)", string=Modifications))
            TheoreticalMonoIsotopicMass, TheoreticalAverageMass = PyMSProcessing.mass_calculator(sequence=Peptide,
                                                                                                 masses_type="MonoIsotopic",
                                                                                                 cysteine_modified=CysteineMod,
                                                                                                 methionine_modified=MethionineMod), PyMSProcessing.mass_calculator(
                sequence=Peptide, masses_type="Average", cysteine_modified=CysteineMod,
                methionine_modified=MethionineMod)


        TheoreticalMZ = PyMSProcessing.mass_to_charge(mass=TheoreticalMonoIsotopicMass, isotope=0, charge_state=Charge)

        # Extracting Ion Chromatogram
        XIC_results = PyMSProcessing.xic(data=Data, mz=ObservedMZ, tolerance=5)

        # Extracting data from XIC
        ## Filtering XIC for only relevant spectra
        XIC_results = XIC_results[XIC_results['Intensity'] > 0]
        ## Option 1: Filtering between a 1 minute window of the original MS2 spectra RT
        #XIC_results = XIC_results[XIC_results['Retention Time'].between(RT_Peptide - 30, RT_Peptide + 30)]
        # Option 2: Just 30 seconds after as the peptide can't be resampled
        XIC_results = XIC_results[XIC_results['Retention Time'].between(RT_Peptide - 5, RT_Peptide + 30)]

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
                            [raw_file, SpectrumID, Peptide, Modifications, PeptideLength, ProteinID, Concentration,
                             PeptideProphetScore, TheoreticalMonoIsotopicMass, TheoreticalAverageMass,
                             Charge, TheoreticalMZ, ObservedMZ, RT, TIC])

                        # Acquiring peak information
                        ## Acquiring the peak information from spectrum
                        array_temp = np.transpose(Spectrum.get_peaks())

                        # Get isotopes
                        Storage_temp = PyMSProcessing.isotope_extractor(peaks=array_temp, storage_list=Storage_temp, mz=ObservedMZ, charge=Charge, tolerance=5)

                        # Add to final storage list
                        Storage.append(Storage_temp)

    # Generate the Pandas dataframe
    Results = pd.DataFrame(Storage, columns=['RawFile',
                                             'Spectrum',
                                             'PeptideSequence',
                                             'Modifications',
                                             'PeptideLength',
                                             'ProteinID',
                                             'Concentration(pmol)',
                                             'PeptideProphetProbability',
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

    Results.to_excel(f"G:/Shared drives/CENSTAT - Bioinformatics/Data/UPS/Article/{raw_file}ObsMZ30Sec.xlsx", index=False)