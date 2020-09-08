# The Paper
This repository includes code for the following paper:

[*Cooper, R.A., Kurkela, K.A., & Ritchey, M. (preprint)*]()

# Abstract

# Code
MRI data was preprocessed with [fMRIPrep](https://fmriprep.org/en/stable/) and denoised with the [CONN toolbox](https://web.conn-toolbox.org/) using fMRIPrep-generated confounds. Functional connectivity analyses were run using Matlab (functions from [cosmomvpa](http://www.cosmomvpa.org/) and [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)) and R Studio.
- Confounds from fMRIPrep preprocessing, including aCompCor, were saved as text files for nuisance regression with the CONN toolbox: See [`preprocessing`](https://github.com/memobc/paper-camcan-pmn/tree/master/preprocessing).
- ROIs within the posterior medial network (PMN) are 100-voxel clusters (2x2x2mm), based on a conjunction of default network subsystems and episodic activity: See [`rois`](https://github.com/memobc/paper-camcan-pmn/tree/master/rois).
- Using denoised data, seed-to-voxel connectivity maps were created for each subject in order to calculate PMN subsystems based on similarity of whole-brain connectivity: See [`analyses`](https://github.com/memobc/paper-camcan-pmn/tree/master/analyses) `wholebrain_seed_to_voxel.m`.
- For intranetwork connectivity analyses in R, the mean time-series of each ROI was extracted from unsmoothed data: See [`analyses`](https://github.com/memobc/paper-camcan-pmn/tree/master/analyses) `extract_roi_timeseries.m`.
- R markdowns run all PMN connectivity analyses using the seed-to-voxel data and intranetwork correlations of ROI time-series — `CamCan-PM-network_[Discovery/Replication].Rmd` — the Group 1 (Discovery) report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Discovery.html) and the Group 2 (Replication) report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Replication.html).
- Finally, we test the relationship between PMN connectivity and episodic memory across all subjects using [intersubject RSA](https://www.sciencedirect.com/science/article/pii/S1053811920303153) — `CamCan-PM-network_Memory.Rmd` — the results report can be viewed [here](http://www.thememolab.org/paper-camcan-pmn/analyses/CamCan-PM-network_Memory.html).

# License
All code is licensed under the [MIT license](https://github.com/memobc/paper-camcanPMN/blob/master/LICENSE).

# Comments?
Please direct any comments or questions to Rose Cooper, rose.cooper at bc.edu. Please feel free to use any of these scripts. Unfortunately I cannot provide support for you to adapt them to your own data. Notice a bug? Please tell me!
