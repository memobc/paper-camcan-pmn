# The Paper
This repository includes code for the following paper:

[**Cooper, R.A., Kurkela, K.A., Davis, S.W., & Ritchey, M. Mapping the organization and dynamics of the posterior medial network during movie watching. *NeuroImage*, doi:10.1016/j.neuroimage.2021.118075**](https://doi.org/10.1016/j.neuroimage.2021.118075)

# Abstract
Brain regions within a posterior medial network (PMN) are characterized by sensitivity to episodic tasks, and they also demonstrate strong functional connectivity as part of the default network. Despite its cohesive structure, delineating the intranetwork organization and functional diversity of the PMN is crucial for understanding its contributions to multidimensional event cognition. Here, we probed functional connectivity of the PMN during movie watching to identify its pattern of connections and subnetwork functions in a split-sample replication of 136 participants. Consistent with prior findings of default network fractionation, we identified distinct PMN subsystems: a Ventral PM subsystem (retrosplenial cortex, parahippocampal cortex, posterior angular gyrus) and a Dorsal PM subsystem (medial prefrontal cortex, hippocampus, precuneus, posterior cingulate cortex, anterior angular gyrus). Ventral and Dorsal PM subsystems were differentiated by functional connectivity with parahippocampal cortex and precuneus and integrated by retrosplenial cortex and posterior cingulate cortex, respectively. Finally, the distinction between PMN subsystems is functionally relevant: whereas both Dorsal and Ventral PM connectivity tracked the movie content, only Ventral PM connections increased in strength at event transitions and appeared sensitive to episodic memory. Overall, these findings reveal PMN functional pathways and the distinct functional roles of intranetwork subsystems during event cognition.

# Code
FMRI data was preprocessed with [fMRIPrep](https://fmriprep.org/en/stable/) and denoised with the [CONN toolbox](https://web.conn-toolbox.org/) using fMRIPrep-generated confounds. Connectivity analyses were run using Matlab (functions from [cosmomvpa](http://www.cosmomvpa.org/) and [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)) and R Studio.
- Confounds from fMRIPrep preprocessing, including aCompCor, were saved as text files for nuisance regression with the CONN toolbox: See [`preprocessing`](https://github.com/memobc/paper-camcan-pmn/tree/master/preprocessing).
- ROIs within the posterior medial network (PMN) are 100-voxel clusters (2x2x2mm), based on a conjunction of the default network and episodic activity: See [`rois`](https://github.com/memobc/paper-camcan-pmn/tree/master/rois).
- Using denoised data, seed-to-voxel functional connectivity maps were created for each subject in order to calculate PMN subsystems based on similarity of whole-brain connectivity: See [`analyses`](https://github.com/memobc/paper-camcan-pmn/tree/master/analyses) `wholebrain_seed_to_voxel.m`.
- For intranetwork connectivity analyses in R, the mean time-series of each ROI was extracted from unsmoothed data: See [`analyses`](https://github.com/memobc/paper-camcan-pmn/tree/master/analyses) `extract_roi_timeseries.m`.
- R markdowns run all PMN connectivity analyses using the seed-to-voxel data and intranetwork correlations of ROI time-series — `CamCan-PM-network_[Discovery/Replication].Rmd` — the Group 1 (Discovery) report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Discovery.html) and the Group 2 (Replication) report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Replication.html).
- Finally, we test the relationship between PMN connectivity and episodic memory across all subjects using [intersubject RSA](https://www.sciencedirect.com/science/article/pii/S1053811920303153) — `CamCan-PM-network_Memory.Rmd` — the results report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Memory.html).

# License
All code is licensed under the [MIT license](https://github.com/memobc/paper-camcanPMN/blob/master/LICENSE).

# Comments?
Please direct any comments or questions to Rose Cooper, rose.cooper at bc.edu. Please feel free to use any of these scripts. Unfortunately I cannot provide support for you to adapt them to your own data. Notice a bug? Please tell me!
