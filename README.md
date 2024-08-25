# Cross-bispectrum decomposition
In this repository you can find code for the decomposition method (`bsfit.m`), the statistical test (`bsfit_stats.m`), simulations (`simulations/`) and data analysis. The method itself has been developed by [Guido Nolte](https://www.uke.de/allgemein/arztprofile-und-wissenschaftlerprofile/wissenschaftlerprofilseite_guido_nolte.html) and [Stefan Haufe](https://www.tu.berlin/uniml/about/head-of-group) and is yet unpublished. 

You can find the pipeline that was used for the simulation experiments in `simulations/main_sim_pac.m`. The main function for the real data analysis is `main.m`. A minimal demo script that shows the decomposition on simulated data without any complex analyses can be found in `main_demo.m`. Running code in this repository requires the installation of [EEGLAB](https://github.com/sccn/eeglab). 

💡 If you have any questions about the code or the project, please reach out to `tien-dung.nguyen@utexas.edu`, I am happy to answer your questions! 

## Problem formulation
The aim of this method is to identify $N$ brain source interactions from given estimates of the cross-bispectrum between $M$ sensors, with $N < M$. Compared to the biPISA approach [Chella et al., 2016](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.93.052420), which decomposes the antisymmetrized sensor cross-bispectrum into a set of pairwise interacting subsystems, this method is not restricted to _pairwise_ source interactions but can be used to evaluate interactions between $n$ sources.

In its most general form, the sensor cross-bispectrum is defined over three channels $i, j, k$ and frequencies $f_1$ and $f_2$ as

$$
B_{ijk}(f_1, f_2) =  \frac{1}{N_e} \ \sum_{e=1}^{N_e} x_{i,e}(f_1) \ x_{j,e}(f_2) \ x_{k,e}^{*}(f_1+f_2),
$$

where $x_{i,e}(f)$ denotes the Fourier transform of the $e^{\text{th}}$ data epoch at frequency $f$ in channel $i$ and $^*$ denotes the complex conjugation. We make the usual assumption that the observed sensor signals $x_i(f)$ result from a linear superposition of the underlying source signals $s_m(f)$, which reads

$$
x_i(f) = \sum_{m=1}^{N} a_{im} \ s_m(f).
$$

The sensor-level cross-bispectrum can therefore expressed as

$$
B_{ijk}(f_1, f_2) = \sum_{l,m,n}^N \ a_{il} \ a_{jm} \ a_{kn} \ D_{lmn}(f_1, f_2),  
$$

with

$$
D_{lmn}(f_1, f_2) = \frac{1}{N_e} \ \sum_{e=1}^{N_e} s_{l,e}(f_1) \ s_{m,e}(f_2) \ s_{n,e}^*(f_1+f_2).  
$$

$D_{lmn}(f_1, f_2)$ can be interpreted as the source-level cross-bispectrum. The aim of the method is to identify subsystems of interacting brain sources from estimates of the sensor cross-bispectrum. This implies finding a set of coefficients $a_{im}$ and a _source_ cross-bispectrum $D_{lmn}(f_1, f_2)$ that together approximate the \textit{sensor} cross-bispectrum. The optimzation problem can therefore be expressed as 


$$
\\{ \boldsymbol{\hat{A}}, \hat{D}_{lmn} \\} =  argmin\_{\\{ \boldsymbol{A}, D\_{lmn} \\}} \ \frac{1}{| B\_{ijk} |} \left|B\_{ijk} - \sum\_{l,m,n}^{N} a\_{il} \ a\_{jm} \ a\_{kn} \ D\_{lmn} \right|
$$

where $\boldsymbol{A}$ is a $M \times N$ matrix that pulls together the individual coefficients $a_{im}$ for channels $M$ and sources $N$, where $N$ is smaller than $M$. The frequency arguments are omitted for ease of reading.
