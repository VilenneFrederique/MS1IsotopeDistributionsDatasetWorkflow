# MS1IsotopeDistributionsDatasetWorkflow
## How to use
Simply download both the PyMSProcessing.py and MS1ProcesssingR.R scripts and place them in the same folder.
Import the PyMSProcessing script and call the necessary functions.

```
import PyMSProcessing
raw_file_test = "A11-12042.mzML"
input_list_test = "Test_list.xlsx"
output_file_test = "Testrun"
PyMSProcessing.excel_MS1_Isotope_Distributions(raw_file=raw_file_test, input_list=input_list_test, output_file=output_file_test)
PyMSProcessing.JSON_MS1_Isotope_Distributions(raw_file=raw_file_test, input_list=input_list_test, output_file=output_file_test)
```

raw_file: An mzML-file with MS1 spectra. The PyMSProcessing script contains a function to extract MS1 spectra from an mzML-file called ms_level_filter, which returns the mzML-file filtered for the MSn spectra specified, which could be either MS1 spectra or MS2 spectra depending on the user-input. 

input_list: The input list with the same format as the Test-list.xlsx in the GitHub Repository

output_file: The name of the output file.



## PyMSProcessing.py
It contains code to extract information from mzML files. It contains a lot of useful functions that may be called separately, such as computing the monoisotopic mass or average mass of different charge states of a peptide given a peptide string and selecting MS1 spectra from an mzML-file.



## MS1ProcesssingR.R
An R-script is called from the PyMSProcessing script to run the BRAIN algorithm.

Dittwald, P., Claesen, J., Burzykowski, T., Valkenborg, D., & Gambin, A. (2013). Brain: a universal tool for high-throughput calculations of the isotopic distribution for mass spectrometry. Analytical Chemistry, 85(4), 1991–1994. https://doi.org/10.1021/ac303439m



## Test_list.xlsx
An example of an input list.



## Citation
When using, please cite:

Vilenne, F., Agten, A., Appeltans, S., Ertaylan, G., & Valkenborg, D. (2023). Establishing a comprehensive workflow for extracting MS1 isotope distributions in LC-MS/MS proteomics. Authorea (Authorea). https://doi.org/10.22541/au.168562296.62840067/v1
