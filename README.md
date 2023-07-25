# CONSISTENT_FDP
This repository contains the code to reproduce all the figures from 
* [False discovery proportion envelopes with consistency.](https://arxiv.org/abs/2306.07819) I. M, G. Blanchard, E. Roquain (2023)*.
 
## Organization of the repository

The code/scripts/ folder contains scripts of the experiments and the method's functions (please consider adapting the file paths if needed).
- For Figure 1, 2, 3, and 4 of Section 5.1 the respective corresponding scripts are topk_simu_vary_m_alpha.R -- topk_sparse_vary_m_alpha.R -- topk_simu_vary_m_alpha.R (with adaptive = 1 in the xp_topk_vary_m_and_alpha.json), and topk_simu_pi0hat.R provided in code/scripts/topk_simu.
- For Figure 5, 6, 7, and 8 of Section 5.2 the corresponding scripts is preordered_simu_vary_m_alpha.R with changes in the parameter .
- For Figure 9
- For Figure 10
<!-- After running the experiments, the raw results are stored in the xp_data/ folder accordingly to the setting of interest. -->
Each time you lauch an experiment a .csv (with date and time detail) file is created in the corresponfind xp_data/ subfolder.
The code/config_files/ folder contains configuration files for the experiments: they allow to set the parameters of the experiments using .json files.
The current parameters are the ones used for the figures provided in the paper.
The code/plots_rmd/ folder contains the Rmd files to reproduce the figures.
At the end of each Rmd, figures are saved to the xp_plot/ folder accordingly to the setting of interest.




<!-- Running the experiments in the main/ folder will provide Figures 6, to 13, 16, and 17 in the figures/ folder and the associated data in the data/ folder.
To launch the experiments type in a terminal
``` 
bash launch_simuxp.sh
bash launch_applixp.sh
``` -->
