## Documentation 
This folder contains the trigger codes and time line used in this experiment. 

### Experiment design 
Participants are presented with a grey screen on which subsequently four items appear clock-wise. 
After all objects appear in their start position, one object moves to another position on a 3x3 grid (see grid_position.png). 
The filler item 1 always moves into the middle (see movement_pattern.xlsx). 
After that a fixation cross appears, followed by a red fixation dot. 
The participants then listens to an audio recording of a male of female speaker explaining the movement they just witnessed (see sentence_pattern.xlsx). 
The first noun phrase is the critical item (see folder lists-and-stimuli for stimuli and experiment lists). 

### Preprocessing 
The folder scripts include all MATLAB preprocessing scripts used in this study. They are divided into a traditional epoch-branch and the regression-based ERP branch (deconvolution branch). 
For this study, we only used the deconvolution branch. 
The traditional epoch-branch is still included for comparison after collecting the whole data set. 

### Statistical analysis 
The folder scripts include all scripts for statistical analyses in MATLAB and R. We used linear mixed models (R) and a cluster-based permutation test (MATLAB).

### Required packages and toolboxes
The required packages for the statistical analysis in R are listed in the corresponding R script. 
The necessary MATLAB toolboxes include EEGLAB, Unfold and FieldTrip. Necessary custom functions are included in the scripts folder. We recommend saving them together as they partially depend on each other. For clarification the scripts are saved to separate folders in the scripts folder. 
