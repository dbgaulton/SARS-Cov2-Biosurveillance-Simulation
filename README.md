# SARS-Cov2-Biosurveillance-Simulation

This repository contains all scripts pertaining to the SARS Cov2 Bisurveillance Simulation Project. Large .fa and zipped .xz and .gz files were omitted from this repository and included in the Large Files release. There are notes in the Large Files Relase that describe what files are necessary for running which .py programs.


## AlternativeMutationalModel
Contains files for the model that propagates mutations based on an entropy threshold. The freq.py file in both the va_delta and va_omicron folders is the python script that generates the mutational model based on the propagation of mutations through a contact network. The output_abridges.py file is the file that defines an abridged contact network from t = 1 to t = 20, which the mutational mode is painted onto.

## hmm
Contains files to generate hmm libraries based on defined genetic cutoffs. The hmm libaries could eventually be used to create an hmm-based mutaitonal model.

### hmm/Global
Contains the files neccessary to build an hmm library based on all COVID sequences globally. Although msaprocessing_global.py should work in theory, the global multiple sequence alignment is too large to run even on a large memory partition. 

### hmm/Virginia/hmm_library
Contains 27 hmms generated based on all Virginia COVID sequences. This libary was generated from the hmm/virginia/msaprocessing_va.py script.

### hmm/Virginia_delta/hmm_library
Contains 27 hmms generated based on Virginia delta COVID sequences. The sequences defined as delta sequences are those that exised while the delta variant was in the majority. That is, from 7/13/2021 to 12/21/2021. This libarary was generated from the hmm/Virginia_delta/msaprocessing_delta.py script.

### hmm/Virginia_omicron/hmm_library
Contains 27 hmms generated based on Virginia omicron COVID sequences. The sequences defined as omicron sequences are those that existed while the omicron variant was in the majority. That is, from 12/21/2021 to 3/7/2022. This library was generated from the hmm/Virginia_omicron/msaprocessing_omicron.py script. 

## network
Contains network.py which generates a visual contact network from t = 0 to 20 (colored by time) using the output_abridges.csv abridged contact network file. 

## trees
Contains newick tree files for all covid sequences globally, in Africa, Asia, Europe, North America, Oceania, and South America. These trees were downloaded from nextstrain.org. These trees can be compared using the robinson-foulds method. 

