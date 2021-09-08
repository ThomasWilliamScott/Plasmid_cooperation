# Plasmid_cooperation

This folder contains supplemental code for the theory component of Dewar et al. (2021), which is described in the main text and in Supp. Info. 4. All code is written for MATLAB R2019a.

The folder contains a script titled "Generate_Data". This script contains code for iterating our recursion equations (given in Supp. Info. 4). All theoretical data presented in Dewar et al. (2021), aside from analytical results, was obtained using this script.

The folder also contains a script titled "Generate_Fig_5", and five data files titled: "Data_s=0.1", "Data_s=0.2", "Data_s=0.3", "Data_s=0.4", "Data_s=0.5". The data files store data obtained by iterating our "Generate_Data" script for different plasmid loss rates (s). The "Generate_Fig_5" script loads these data files and uses them to plot Fig. 5 in the main text.
