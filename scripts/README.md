Here you will find all relevant MATLAB and R scripts for the preprocessing and analysis of the collected EEG data. 

The MATLAB scripts include a general preprocessing pipeline which processes each participant separately and creates files accordingly for each processing steps, a preprocessing script for traditional epoching and reject as well as a preprocessing script for regression-based EPR estimation to better account for overlaps in the neural signals. 

Following these two scripts, there are scripts for exporting GAMM data to use in R and cluster-based permutation tests. One preprocessing branch including export and cluster-based permutation should be run fully before switching to the other preprocessing branch. The MATLAB scripts include all custom functions used for preprocessing and analysis. 

The R scripts include an rmd.-file for running a sensitivity analysis to determine the optimal sample size as well as analyses scripts, including a generalised linear model and generalised additive mixed model approach. 

Please notice: all Matlab files should be stored together as some use the same function-files that are currently stored within different folders. 
Please prepare your Analysis folder the following way: 
Analysis
├── Raw data
│   ├── VP1_raw_data.vdhr
│   ├── VP1_raw_data.vmrk
│   ├── VP1_raw_data.txt
│   ├── VP1_raw_data.eeg
|    ......
├── matlabscript1.m
├── matlabscript1.m
...
├── Preprocessed
│   ├── ...
├── Unfold_Preprocessed
│   ├── ...
├── channel locations 
│   ├── CMA-26-X1.bvef

Matlab will find and assign old and new files to the correct folders. For the Deconvolution it might prove necessary to move the Px_ica_preprocessed.set files to the Unfold_Preprocessed folder or adjust the script to that Matlab can find the correct files. 
