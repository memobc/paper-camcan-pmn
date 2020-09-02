# The Paper
This repository includes scripts and data for the following paper:

# Abstract

# Analysis Code
MRI data was preprocessed with [fMRIPrep](https://fmriprep.org/en/stable/) and denoised with the [CONN toolbox](https://web.conn-toolbox.org/) using fMRIPrep-generated confounds. Functional connectivity analyses were run using Matlab (functions from [cosmomvpa](http://www.cosmomvpa.org/) and [SPM](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)) and R Studio.
- Confounds from the fMRIPrep preprocessing pipeline, including aCompCor, were saved as text files for nuisance regression in the CONN toolbox: See [preprocessing](https://github.com/memobc/paper-camcanPMN/tree/master/preprocessing).
- ROIs within the posterior medial network (PMN) are 100-voxel clusters (2x2x2mm), based on a conjunction of default network subsystems and episodic activity: See [rois](https://github.com/memobc/paper-camcanPMN/tree/master/rois).
- Using denoised data, seed-to-voxel connectivity maps are created for each subject in order to calculate PMN subsystems based on similarity of whole-brain connectivity: See [analyses](https://github.com/memobc/paper-camcanPMN/tree/master/analyses) `wholebrain_seed_to_voxel.m`.
- 

# License
All code is licensed under the [MIT license](https://github.com/memobc/paper-camcanPMN/blob/master/LICENSE).

# Comments?
Please direct any comments or questions to Rose Cooper, rose.cooper at bc.edu. Please feel free to use any of these scripts. Unfortunately I cannot provide support for you to adapt them to your own data. Notice a bug? Please tell me!
