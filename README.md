# Cross-bispectrum decomposition
In this repository you can find code for the decomposition method (`bsfit.m`), the statistical test (`bsfit_stats.m`), simulations (`simulations/`) and data analysis. The method itself has been developed by [Guido Nolte](https://www.uke.de/allgemein/arztprofile-und-wissenschaftlerprofile/wissenschaftlerprofilseite_guido_nolte.html) and [Stefan Haufe](https://www.tu.berlin/uniml/about/head-of-group) and is yet unpublished. 

Running code in this repository requires the installation of [EEGLAB](https://github.com/sccn/eeglab) and its plugin [ROIconnect](https://github.com/sccn/roiconnect?tab=readme-ov-file). Please make sure to have the most up to date version of ROIconnect installed.

## Problem formulation
The aim of this method is to identify $n$ brain source interactions from given estimates of the cross-bispectrum between $k$ sensors, with $n < k$. Compared to the biPISA approach [Chella et al., 2016](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.052420), which decomposes the antisymmetrized sensor cross-bispectrum into a set of pairwise interacting subsystems, this method is not restricted to _pairwise_ source interactions but can be used to evaluate interactions between $n$ sources.

Given the 
