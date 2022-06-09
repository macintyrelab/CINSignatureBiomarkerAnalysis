## Copy-number signatures as biomarkers for drug response prediction and drug target identification

This repository contains the code for reproducing Figure 3 of Drews et al. (2022). The central hub for the study can be found in https://github.com/markowetzlab/Drews2022_CIN_Compendium

Run scripts to reproduce the analysis focused on testing CINsignatures as biomarkers for drug response and discovery of novel drug targets

#### Summary
1) Correlation analyses between CIN signature activities, CRISPR knock-out screens, RNAi knock-down screens, and PRISM drug responses (run trifecta_analysis.Rmd). Ouputs of this script are included in the data folder. They are needed to run the DrugScreening.R script
3) Response biomarkers: combination of significant lists: CRISPR+drug, RNAi+drugs, and CRISPR+RNAi+drugs --> look for drug's MOA 
4) Novel drug targets: combination of CRISPR+RNAi without known drug --> selection of druggable targets according to canSAR
5) Alluvial plots

## Contact

If you experience any issues or have questions about the code, please open a Github issue with a minimum reproducible example. For questions around collaborations or sensitive patient data, please contact us directly at Florian Markowetz <Florian.Markowetz@cruk.cam.ac.uk> and Geoff Macintyre <gmacintyre@cnio.es>.

## Licence
The contents of this repository are copyright (c) 2022, University of Cambridge and Spanish National Cancer Research Centre (CNIO).

The contents of this repository are published and distributed under the GAP Available Source License v1.0 (ASL). 

The contents of this repository are distributed in the hope that it will be useful for non-commercial academic research, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ASL for more details. 

The methods implemented in the code are the subject of pending patent application GB 2114203.9.

Any commercial use of this code is prohibited.
