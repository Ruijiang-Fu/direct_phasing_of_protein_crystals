# direct_phasing_of_coiled_coil_protein_crystals
Software requirement:
1. Intel MKL libraries should be installed. Intel MKL is part of Intel One API which can be downloaded free 
from https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html
2. When compile the source code, the following link should be used 
    -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm -O

Data source:
The observed data in cif and the posted structure in pdb are downloaded from the protein data bank at www.rcsb.org

Notes for the source code in each directory:
1. catch_observed_reflections_from_cif_file is used to catch reflections from the downloaded cif file.
2. arrange_unique_reflections is used to sort the reflections by resolution.
3. calculate_boundary_by_FFT is used to calculate a protein boundary from downloaded pdb.
4. calculate_standard_histogram is used to calculate a reference histogram.
5. calculate_correct_phases_for_different_origin_choices is used to calculate the correct phases for different origin choices.
6. directPhasing is used to retrieve the lost phases. The correct phases in that directory is used for computing the mean phase error.

Readers can try running the code in direct_phasing_of_3ILZ or direct_phasing_of_1Y5Y. When apply the cource code on new diffraction data, it requires some experience to acheive a solution. Each run in directPhasing can be compiled and run on a desktop, but a multi-core computing server is preferred where tens of runs can be compiled and run at the same time.
We are trying to improve the source code and will update new versions. If you have unphased diffraction data, 
you could reach us via email. It's our pleasure to help you.
