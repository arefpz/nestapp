# nestapp v2.0
The NeuroStimulation App (nestapp) utilizes EEGLAB and TESA functions to clean the TMS-EEG data automatically and plot the TMS evoked potentials (TEPs) for selected or already cleaned data.
The app now has two tabs: Cleaning and Visualizing.
# Introduction
The nestapp has been developed using MATLAB App Designer (R2023b & R2018a). The idea is to use the currently available EEGLAB (v2025.0.0), TESA (v1.1.1), and FastICA (2.5) functions for cleaning and analyzing TMS-EEG data.
Note that FastICA should be added to your MATLAB path.
# How to use
In the Cleaning tab:
1- After downloading the app, run nestapp.m. 

2- On the left panel, you will see a list of steps that you can choose to apply to your data. By clicking on each step, extra information including the function name, will be available in the 'Info' section. Please read the extra information about each function, on EEGLAB, TESA, and FastICA websites:

EEGLAB: https://eeglab.org/
TESA: https://nigelrogasch.gitbook.io/tesa-user-manual/
FastICA: http://research.ics.aalto.fi/ica/fastica/index.shtml

3- You can use the 'Add', 'Remove', 'Move Up' and 'Move Down' buttons to create your desired cleaning pipeline in the 'Selected Steps'. Each of the steps may have parameters that you wish to modify. You can change the value for appropriate properties in the panel below the 'Parameter(s)' section by clicking on the value for appropriate properties. The 'Default Value' button will reset the values of the selected step into the default values.

4- You can save the created pipeline for future uses, by clicking 'Save pipeline' and giving it the desired name. The pipeline will be saved as a MATLAB Mat file. The 'Load Pipeline' will load the pipeline into the 'Selected Steps' section.

5- In the right panel, if EEGLAB is not in your MATLAB path, you can select the folder containing the EEGLAB package.

6- In the 'Select data to perform Analysis' section, choose the files ( in *.vhdr, *.set, *.cnt, *.cdt format) you want to apply the analysis. the number of selected files and the index of the current file being processed will appear below.

7- Pressing 'Re/Start Steps' will erase all selected steps and changes made to the parameters.

8- By pressing the 'Run Analysis' the cleaning and analyzing process will be started.

NOTE: In case your data doesn't have a channel location and you would like to use 'Load Channel Location', please make sure that this function is performing right. Extra features will be added in the future.

In the Visualizing tab:
Now, you can plot the TEPs for selected data. To plot the TEPs, you need to select the Region of Interest (ROI) electrodes by clicking on them. The selected electrodes will be highlighted. Using the buttons, you can either plot the selected TEP in a new figure or add it to the existing plot. You can also export the TEP for the selected files with a given name to the MATLAB Workspace for further analysis. Also, the topoplot button will allow you to plot the overall activity of the brain around the selected time point and averaged over a given time window. You can also plot the EEG data by selecting the data from the drop-down list.

# Code Contributors/Developers
Aref Pariz: The developed code (v1.0 2023) was a part of my job at the Royal Institute for Mental Health, in Dr. Sara Tremblay's lab (NESTLAB: https://www.nest-lab.ca/) and Dr. Jeremie Lefebvre's Lab at the University of Ottawa.

