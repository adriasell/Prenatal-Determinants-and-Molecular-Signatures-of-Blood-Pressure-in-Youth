# From Environment to Biology: Prenatal Determinants and Molecular Signatures of Blood Pressure in Youth
This repository contains the complete analytical pipeline and custom code to reproduce the findings presented in the study "From Environment to Biology: Prenatal Determinants and Molecular Signatures of Blood Pressure in Youth".

## Overview
The study explores the interplay between prenatal environmental exposures, molecular signatures, and adolescent cardiovascular health. We integrate multiple omics layers to identify biomarkers of blood pressure (BP) and their upstream environmental drivers.
## System Requirements
* Memory: Minimum 16GB RAM (32GB+ recommended for large multi-omic matrices).
* R version: 4.4.1
## Data availability
The data supporting the findings from this study are available within the manuscript and its supplementary information. Due to the HELIX data policy and data use agreement, human subjects’ data used in this project cannot be freely shared. Researchers external to the HELIX Consortium who have an interest in using data from this project for reproducibility or in using data held in the HELIX data warehouse for research purposes can apply for access to data for a specific manuscript at the time. Interested researchers should fill in the application found at https://athleteproject.eu/helix-cohort/ and send it to helixdata@isglobal.org. The applications are received by the HELIX Coordinator, and are processed and approved by the HELIX Project Executive Committee. The decision to accept or reject a proposal is taken by the HELIX Project Executive Committee, and is based largely on potential overlap with other HELIX-related work, the adequacy of data protection plans, and the adequacy of authorship and acknowledgement plans. Further details on the content of the data warehouse (data catalogue) including those data used for the present study and procedures for external access are described on the project website https://athleteproject.eu/helix-cohort/.
## Repository Structure
The repository is organized to distinguish between the core analytical pipeline, custom utility functions, and the specific software environment used:
* Codi: Contains the primary R scripts used to execute the multi-omics integration and environmental driver analysis.
* Packages: Contains a modified version of the RGCCA R package. This version allows to set sparsity to 0, which essential for the specific requirements of this study.
* Scripts: Contains auxiliary R functions and custom utility scripts called by the main analysis files.
