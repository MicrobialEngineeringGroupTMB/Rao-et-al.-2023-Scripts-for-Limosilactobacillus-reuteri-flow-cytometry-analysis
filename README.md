# Rao et al. (2023) -  Scripts for *Limosilactobacillus reuteri* flow cytometry analysis

## Background

The scripts stored in this repository contains the code used for the processing of flow cytometry data for the Rao et al (2023) publication titled *Non-inhibitory levels of oxygen during cultivation increase freeze-drying stress tolerance in Limosilactobacillus reuteri DSM 17938*. The scripts were developed using data from the bacterium Limosilactobacillus reuteri, which was stained with SYBRgreen and Propidium iodide and captured with a BD Accuri C6+ flow cytometer. 

The repository contains three Matlab scripts, all of which are designed to process all .fcs files in a user-specified directory:

* a wrapper script (*FCSData_to_MData*) to call on a fcs reader (Balkay, 2023) to load the flow cytometry data (FCS3.1 format) into Matlab

* a script for fixed gating (*FCM_analysis_live_dead_damaged_fixed_gates_2023*)

* a script for automated gating using k-means clustering (*FCM_analysis_bacteria_kmeans_2023*)

Notebooks with examples data are available for the <a href="https://microbialengineeringgrouptmb.github.io/Rao-et-al.-2023-Scripts-for-Limosilactobacillus-reuteri-flow-cytometry-analysis/fixed_gates_live_script.html">fixed gating script</a> and for the k-means gating script.

If there is anything in this repository that you found useful for your own work, we would be grateful if you would consider citing the Rao et al (2023) paper.

## Dependencies

* Matlab (the initial commit was tested with version 2022a) 

* fca_readfcs (version used for the initial commit was 2020.06.22)
https://se.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs


* the following Matlab toolboxes: Statistics and Machine Learning, Curve Fitting, Fixed-point Designer

(NB! There are multiple functions within the official Matlab toolboxed named *nearest*, and the one used in these scripts is the one that is included in the Fixed-point Designer toolbox)

## References

 Balkay, L. (2023). fca_readfcs ([https://www.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs](https://www.mathworks.com/matlabcentral/fileexchange/9608-fca_readfcs)), MATLAB Central File Exchange. Retrieved January 20, 2023. 

 Rao, N. S., Lundberg, L. E., Tomasson, J., Tullberg, C., Brink, D.P., Lundberg, L., Palmkron, S. B., van Niel, E. W. J., Håkansson, S., Carlquist, M.(2023). Non-inhibitory levels of oxygen during cultivation increase freeze-drying stress tolerance in *Limosilactobacillus reuteri* DSM 17938. (Submitted)
