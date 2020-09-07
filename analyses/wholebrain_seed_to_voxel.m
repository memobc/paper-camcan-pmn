%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using pre-defined ROIs, this script runs wholebrain seed-to-voxel
% connectivity analyses.
%
% A connectivity map of fisher-z transformed r values is saved out 
% for each seed (ROI).
%
% Note - seed data is unsmoothed and voxel data is smoothed (6mm).
% Denosied data used here is already masked with an MNI brainmask.
%
% Rose Cooper - Feb 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
warning('off','all');

b.scriptdir = '/gsfs0/data/ritcheym/data/fmri/CamCan/pm-network/analysis/scripts/';
addpath(genpath(b.scriptdir));

% toolboxes
b.spmDir = '/data/cooperrn/Documents/fmri-core/'; %spm/conn/cosmo
addpath(genpath(b.spmDir));


%%% define task to analyse
b.task = 'movie';
% group to analyse? Large data-set -- used split-sample replication
b.sample = 'Replication'; %'Discovery';
%%%%%%%%%%%%%%%%%%%%%%%%%%


base_dir = '/gsfs0/data/ritcheym/data/fmri/CamCan/pm-network/';

%where are my ROI files?:
b.ROIdir   = [base_dir 'analysis/rois/pm_rois/'];
roi_names  = table2cell(readtable([b.ROIdir 'PM_voxel_clusters_IDs.csv']));
% ^^ col 1 = name, col 2 = ID in nifti
roi_labels = cosmo_fmri_dataset([b.ROIdir 'rPM_voxel_clusters.nii']);
roi_labels = roi_labels.samples;


%where is the functional data?:
b.dataDir = [base_dir 'data/derivs/denoised/'];
%save .nii files to?:
b.outDir  = [base_dir 'analysis/seed-to-voxel/' b.sample '/'];
if ~exist(b.outDir,'dir')
    mkdir(b.outDir);
end


%define subjects to analyse *in this group only*
subjs = table2cell(readtable([base_dir 'data/valid_camcan_network-analysis_groups.csv'], 'ReadVariableNames', 1));
subjs = subjs(contains(subjs(:,4),b.sample),1);


%% First-level connectivity maps

fprintf('\nCalculating seed to voxel (whole brain) connectivity for %s sample... \n\n',b.sample);
for s = 1:length(subjs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curSubj = subjs{s};
    
    fprintf('\t%s\n',curSubj);
    b.subDir = [b.outDir curSubj '/'];
    mkdir(b.subDir);
    

    % use cosmo to grab data: -------------------------------------
    % load ROI data (unsmoothed)
    myFile = [b.dataDir curSubj '_task-' b.task '_MNI_denoised.nii'];
    template = spm_vol([myFile ',1']); % for output space
    node_data = cosmo_fmri_dataset(myFile);
    % set all zeros to NaNs (should be NaNs already)
    node_data.samples(node_data.samples == 0) = NaN;
 
    % load brain data (smoothed)
    myFile = [b.dataDir 'sm_' curSubj '_task-' b.task '_MNI_denoised.nii'];
    brain_data = cosmo_fmri_dataset(myFile);
    % set all zeros to NaNs (should be NaNs already)
    brain_data.samples(brain_data.samples == 0) = NaN;
    

    % run seed to voxel for each ROI
    for r = 1:size(roi_names,1)
        
        fprintf('\t\tworking on %s, %s\n', b.task, roi_names{r,1});
        
        %get roi data (nTR x nVox):
        roi_data = node_data.samples(:,roi_labels == roi_names{r,2});
        if sum(isnan(roi_data(:))) > 0
            error('NaNs present in data');
        end
        tr_means = mean(roi_data,2); % mean across cols (voxels) per row (trial/tr)
        
        % run seed-to-voxel correlations (pearson's r, converted to z). 
        % Note values will be NaN if outside brainmask in whole brain data
        conn_values = atanh(corr(tr_means,brain_data.samples));       

        % save .nii for this roi and task
        v = template;
        v.fname = [b.subDir curSubj '_task-' b.task '_' roi_names{r,1} '_Z_seedtovoxel.nii'];
        conn_values = reshape(conn_values, v.dim);
        spm_write_vol(v, conn_values); 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end %end of loop through subjects



%% Group-averaged connectivity maps
 
% list all seed-to-voxel first level files (1 per subj, per ROI)
map_files = struct2cell(dir([b.outDir '**/']));

b.grpOut = [b.outDir 'group-average/'];
if ~exist(b.grpOut, 'dir')
    mkdir(b.grpOut);
end

roi_patterns = [];
for r = 1:length(roi_names)
    
    fprintf('\n\nCalculating group-average for %s, %s\n', b.task,roi_names{r,1});
    
    % grab subject files for this task and roi
    current_files = map_files(1:2, contains(map_files(1,:),b.task));
    current_files = current_files(1:2,contains(current_files(1,:),roi_names{r,1}));
    current_files = current_files(1:2,contains(current_files(1,:),'sub-'));
    
    fprintf('\tAveraging across %d subjects\n',size(current_files,2));
    
    
    % now load data for each subject
    template = spm_vol([current_files{2,1},'/',current_files{1,1}]);
    maps = nan([template.dim length(current_files)]);        
    for s = 1:size(current_files,2)
        maps(:,:,:,s) = spm_read_vols(spm_vol([current_files{2,s},'/',current_files{1,s}]));
    end
    
    % take the average across subjects, convert back to r, and save
    map_average = tanh(mean(maps,4));
    roi_patterns(:,r) = map_average(~isnan(map_average)); % correlations to non-nan voxels
    
    % save .nii for this roi and task
    template.fname = [b.grpOut 'task-' b.task '_' roi_names{r,1} '_R_seedtovoxel.nii'];
    spm_write_vol(template, map_average);
end

% save out group-averaged csv with voxel connectivity values per roi
output = array2table(roi_patterns);
output.Properties.VariableNames = roi_names(:,1);
outName = [b.grpOut 'group_PM_wholebrain_connectivity.csv'];
writetable(output, outName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%