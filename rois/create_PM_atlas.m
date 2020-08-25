%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script defines posterior medial network ROIs according to
% the following criteria (as defined in Ritchey & Cooper, 2020, TICS):
%
% 1. Regions from Schaefer et al. (2018), that fall within Default A and C
% subnetworks, and posterior hippocampus (Ritchey et al., 2015), all 
% left hemisphere.
%
% 2. Conjunction between the DMN regions (above) and voxels showing 
% "episodic" sensitivity from a Neurosynth meta-analysis.
%
% From these masks, 100-voxel clusters are selected (2x2x2 space) 
% expanding from peaks of episodic activity.
%
%%%%%%%%%%%%%%%%%%%
% Rose Cooper, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
warning('off','all')

target_size = 100; %voxel cluster size for pm regions


% --------------------------------------------------------------------- %
% step 1 --> filter schaefer 2018 parcels by those in (left) MTL and core
% default systems (labelled as A and C subnetworks)
vol = spm_vol('Schaefer2018_400Parcels_17Networks_order_FSLMNI152_2mm.nii');
parcels     = spm_read_vols(vol);
DMN         = table2cell(readtable('Schaefer2018_400Parcels_MTL-Core.csv','ReadVariableNames',0));
[~,parcels_DMN] = ismember(parcels,cell2mat(DMN(:,1))); %returns rois, now with rescaled ids
% --------------------------------------------------------------------- %


% --------------------------------------------------------------------- %
% step 2 --> add post. hipp regions to schaefer parcels (already resliced to
% Schaefer et al. space)
hipp  = spm_read_vols(spm_vol('HIPP_BODY_L_mask.nii')) + ...
        spm_read_vols(spm_vol('HIPP_TAIL_L_mask.nii'));

top_idx = max(unique(parcels_DMN(:)));
parcels_DMN(logical(hipp)) = top_idx + 1;
% --------------------------------------------------------------------- %


% --------------------------------------------------------------------- %
% step 3 --> re-label regions in the parcellation according to where they
% are (in line with Schaefer et al. descriptive labels):
% accounts for network division and broad region:
regions = unique(DMN(:,3)); % unique combo of network + area

% first add hipp, then cycle through cortical regions to number
parcels_PM = zeros(size(parcels_DMN));
parcels_PM(ismember(parcels_DMN,top_idx + 1)) = 1; %HIPP
for r = 1:length(regions)
    parcels_PM(ismember(parcels_DMN,cell2mat(DMN(contains(DMN(:,3),regions{r}),4)))) = r+1;
end
regions = [{'Hipp'};regions];

% save full PM mask, for reference, to refine by episodic activity (below)
vol.fname = 'PM_mask.nii';
spm_write_vol(vol,parcels_PM);
% --------------------------------------------------------------------- %


% --------------------------------------------------------------------- %
% step 4 --> get "episodic" voxels from Neurosynth (already resliced)
[data, xyz] = spm_read_vols(spm_vol('Neurosynth_episodic_association-test_z_FDR_0.01.nii'));

% a) mask with Schaefer DMN subsystem (+ hipp) regions:
data(parcels_PM == 0) = 0;

% b) cluster-threshold (my function) to remove very small isolated peaks
% min X voxels, adjoining faces (6), min value
data = cluster_threshold(data, 100, 6, 0);

% save DMN subnetworks (+ hipp) and episodic conjunction mask for reference
pm_mask = parcels_PM;
pm_mask(data == 0) = 0;
vol.fname = 'PM_episodic_mask.nii';
spm_write_vol(vol,pm_mask);

% c) find local peaks in episodic sensitivity:
% use all dimensions, min height(z), min voxel distance
% only look for multiple peaks in MPC due to lack of differentiation
% between Prec and PCC
coords = [];
region_ids = unique(parcels_PM(data > 0))'; %in case any lost through above thresholding
for r = region_ids
   mask = data;
   mask(parcels_PM ~= r) = 0;
   if r == 2, th = 6; else, th = max(mask(mask > 0)); end

   peaks = findpeaksn(mask, true(1, ndims(mask)), th, 12);
   v = xyz(:,peaks(:));
   coords = [coords, [v;repmat(r,1,size(v,2))]];
end
% --------------------------------------------------------------------- %


% --------------------------------------------------------------------- %
% step 5 --> remove any peaks too close across regions
new_coords = coords;
for c = 1:(size(coords,2)-1)
    idx = (c+1):size(coords,2);
    dist = sqrt(sum((coords(1:3,c) - coords(1:3,idx)) .^ 2));
    if min(dist) < 12
        new_coords(:,c) = nan(4,1);
    end
end
new_coords(:,isnan(new_coords(1,:))) = [];
% --------------------------------------------------------------------- %


% --------------------------------------------------------------------- %
% step 6 --> create voxel clusters by expanding, within regional mask, from
% each peak
pm_clusters = zeros(size(parcels_PM));
cluster_coords = [];
for c = 1:size(new_coords,2)

    %define vector for pm_mask (episdic voxels in this region)
    %based on the coordinate's region
    mask_vector = logical(pm_mask(:) == new_coords(4,c));

    while 1
        % now allocate voxels values based on their distance from the peak:
        values = abs(xyz - new_coords(1:3,c));
        values = sum(values,1);
        values(2,:) = 1:length(values); %add indexing for re-sorting

        %set voxels out of mask to be invalid
        values(:,~mask_vector) = [];

        % grab X closest voxels to peak *within mask of that region*
        values = sortrows(values',1);
        new_mask = zeros(size(pm_mask(:)));
        new_mask(values(1:target_size,2)) = 1; %mark closest X voxels as 1 based on indexed position.

        % reshape to 3D matrix
        ROImap = reshape(new_mask,size(pm_mask));

        %%% test to make sure that we have a contiguous cluster, with at
        %%% least faces touching
        clusters = bwconncomp(ROImap,6);
        if length(clusters.PixelIdxList) ==1
            pm_clusters(clusters.PixelIdxList{1}) = c;
            cluster_coords(:,c) = mean(xyz(:,logical(new_mask)),2);
            break;
        else %if not a contiguous cluster, exclude smallest cluster and repeat
            lengths = cellfun('length', clusters.PixelIdxList);
            min_length = min(find(lengths == min(lengths)));
            mask_vector(clusters.PixelIdxList{min_length}) = 0;
        end
    end
    
    % now save as an individual binary mask:
    this_roi = zeros(size(pm_clusters));
    this_roi(pm_clusters == c) = 1;
    vol.fname = ['PM_voxel_cluster_' num2str(c) '.nii'];
    spm_write_vol(vol,this_roi);
end
disp(cluster_coords);
fprintf('N clusters = %d\n',size(cluster_coords,2));
% --------------------------------------------------------------------- %

% save pm_clusters file for all regions
vol.fname = 'PM_voxel_clusters.nii';
spm_write_vol(vol,pm_clusters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%