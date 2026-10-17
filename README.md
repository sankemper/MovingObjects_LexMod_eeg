# MovingObjects_LexMod EEG Study
Analysis scripts and reproducibility materials for a Master's thesis on gender informativity in referential processing using EEG

## Overview
This repository contains the analysis scripts and documentation
for a Master's thesis investigating listeners' sensitivity to
gender informativity during spoken German language comprehension.

The study examines the effects of gender informativity,
prenominal modification, and mention on neural responses
to critical nouns.

## Repository structure
- `scripts/`: MATLAB and R analysis scripts.
- `stimuli-and-lists/`: Approved experimental materials and lists.
- `documentation/`: Experimental design, analysis workflow,
  and software requirements.

## Analysis workflow
The analysis pipeline includes:

1. EEG preprocessing and ICA component rejection.
2. Trial relabelling and event timing checks.
3. Deconvolution analysis using unfold.
4. Epoching and artefact rejection.
5. ERP analysis, including predefined regions of interest.
6. Statistical analyses using GAMMs and cluster-based
   permutation tests, where applicable.

See the documentation for the exact order of execution
and the status of each analysis stage.

## Software and dependencies
The analysis uses MATLAB and R, together with relevant
toolboxes and packages.

## Data availability
Participant-level EEG recordings and other restricted
research data are not included in this repository.

## Reproducibility
The repository documents the analysis procedures and
software dependencies. Reproducing the results may require
access to the original research data, which are not
distributed here.
