%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracts functional time series (unsmoothed) from a list of ROIs and list
% of subjects
%
% Rose Cooper Feb 2020
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
nTR    = 193;
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
b.outDir  = [base_dir 'analysis/mean-timeseries/'];
if ~exist(b.outDir,'dir')
    mkdir(b.outDir);
end


%define subjects to analyse *in this group only*
subjs = table2cell(readtable([base_dir 'data/valid_camcan_network-analysis_groups.csv'], 'ReadVariableNames', 1));
subjs = subjs(contains(subjs(:,4),b.sample),1);


% to store time-series in long format 
timeseries = num2cell(zeros([(nTR*size(roi_names,1)*length(subjs))+1,4]));
timeseries(1,:) = {'Subject','Time','Node','Value'};

rowN = 1; %to count rows for timeseries matrix
for s = 1:length(subjs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    curSubj = subjs{s};
    
    fprintf('\nExtracting mean time-series data for %s... \n',curSubj);
    
    
    % use cosmo to grab data: ------------------------------------- 
    % Load data into cosmo (unsmoothed file, no 'sm_' prefix)
    myFile = [b.dataDir curSubj '_task-' b.task '_MNI_denoised.nii'];
    data = cosmo_fmri_dataset(myFile);
    % should already be NaNs outside brainmask but, in case not, set all zeros to NaNs
    data.samples(data.samples == 0) = NaN;
    
    % loop over rois
    for r = 1:size(roi_names,1)    
        % mask with roi (nTR x nVox)
        roi_data = data.samples(:,roi_labels == roi_names{r,2});
        
        if sum(isnan(roi_data(:))) > 0 %any missing voxels in mask?
            error('NaNs present in ROI data');
        end
        tr_means = mean(roi_data,2); % mean across cols (voxels) per row (trial/tr)
        
        % add to output:
        curRows = (rowN+1):(rowN+length(tr_means));
        
        timeseries(curRows,1) = {curSubj};
        timeseries(curRows,2) = num2cell(1:length(tr_means));
        timeseries(curRows,3) = roi_names(r,1);
        timeseries(curRows,4) = num2cell(tr_means);
        
        rowN = rowN + length(tr_means);  
        
        fprintf('\t%s, nvox = %d\n',roi_names{r,1},size(roi_data,2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------
fprintf('Done! Saving timeseries to csv file ... \n\n');

% save CSV file:
resultsTable = array2table(timeseries(2:end,:));
resultsTable.Properties.VariableNames = timeseries(1,:);
outName = [b.outDir 'PM_node_timeseries_' b.sample '.csv'];
writetable(resultsTable, outName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%