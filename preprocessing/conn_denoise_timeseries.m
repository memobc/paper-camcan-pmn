function conn_denoise_timeseries(j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function denoises functional time series by running CONN's confound
% regression pipeline
%
% conn saves the denoised voxel-level 4D data in:
% /results/preprocessing/niftiDATA_Subject*_Condition*.nii
% This file is moved and renamed in ~data/derivs/denoised/
%
% This function runs denoising for a single subject (stored in their own
% project) for ease of parallel processing. Alternatively, subjects could
% be combined into 1 project. Note that conn also has internal options to
% manage parallel processing of subjects, but it will try to unzip all 
% subject niftis iteratively beforehand, which takes a long time. 
%
% help(conn_batch) for options
%
%
% Rose Cooper - Feb 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning('off','all');
b.scriptdir = '/gsfs0/data/ritcheym/data/fmri/CamCan/pm-network/analysis/scripts/';
addpath(genpath(b.scriptdir));


%%% define task to denoise
b.task      = 'movie';
b.nsessions = 1;  % runs
%%%%%%%%%%%%%%%%%%%%%%%%%%


% SPM/CONN info
b.spmDir = '/data/cooperrn/Documents/fmri-core/'; %spm and conn
% Note that my spm_defaults' (saved as *spm_my_defaults*) implicit threshold
% is set to -inf - so we just use (and need) an explicit mask specified.
addpath(genpath(b.spmDir));


base_dir = '/gsfs0/data/ritcheym/data/fmri/CamCan/pm-network/data/';

%where is preprocessed data?
b.derivDir = [base_dir 'derivs/fmriprep/'];
%where are my confound regressors?
b.covDir   = [base_dir 'derivs/confounds/'];
%where to put denoised data?
b.analysisDir  = [base_dir 'derivs/denoised/'];
if ~exist(b.analysisDir,'dir')
    mkdir(b.analysisDir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear batch


%% 1. SETUP SUBJECT(S)

batch.Setup.isnew     = 1;    %0 if want to update existing project
batch.Setup.RT        = 2.47; %TR (seconds)
batch.Setup.done      = 1;    %0 to specify files files but not run step
batch.Setup.overwrite = 1;    %0 if updating existing project and want to keep previous files


% load in file with *valid* subject IDs (to include in analyses):
b.subjs = table2cell(readtable([base_dir 'valid_subjects_camcan-pmnetwork.csv'], 'ReadVariableNames', 0));
b.subjs = b.subjs(j); % running in parallel as separate jobs, so select current subject
% ^^ comment this last line out for running multiple subjects together


batch.Setup.nsubjects = length(b.subjs);
batch.filename=fullfile(b.analysisDir,[b.subjs{1} '_task-' b.task '.mat']); % New conn_*.mat project name
% ^^ here, I'm naming by subject, but could just name by task for subjects processed together


% now only run if we haven't already denoised data for this subject
if ~exist(batch.filename,'file')  
    %% 2. GET FUNCTIONALS
    
    fprintf('\nGetting functional scans...\n');
    
    batch.Setup.functionals = repmat({{}},[1,1]); % initialize main functional volumes for each subject/session
    
    % grab all niis in preprocessed functional folder, as only valid runs were copied there
    % a) main 4D functionals :
    for nses = 1:b.nsessions
        for nsub = 1:length(b.subjs)
            func_regexp = ['^' b.subjs{nsub} '.*' b.task '.*MNI.*preproc.*.nii.gz']; %fmriprep output (unsmoothed)
            scanRuns = cellstr(spm_select('FPList', [b.derivDir b.subjs{nsub} '/func/'], func_regexp)); % full path
            if isempty(scanRuns{1})
                error('no functional files found!');
            end
            batch.Setup.functionals{nsub}{nses} = scanRuns{1};
            
            % unzip functionals - for now, conn requires .nii
            gunzip(scanRuns);
        end
    end
    
    
    %% 3. GET STRUCTURALS
    
    %technically I don't think conn needs these files because I ran
    %aCompCor through fmriprep, but adding them just in case
    
    fprintf('\nGetting structural scans...\n');
    
    batch.Setup.structurals = repmat({{}},[1,1]); %initialize variable for structural scan names
    
    % find MNI space T1 preprocessed image from fmriprep
    for nsub = 1:length(b.subjs)
        % T1
        strc_regexp = ['^' b.subjs{nsub} '.*MNI.*preproc_T1w.nii.gz'];
        strcFile = cellstr(spm_select('FPList', [b.derivDir b.subjs{nsub} '/anat/'], strc_regexp));
        if isempty(strcFile{1})
            error('no T1 file found!');
        end
        batch.Setup.structurals{nsub} = strcFile{1};
        
        % unzip T1 - for now, conn requires .nii
        gunzip(strcFile);
        
        
        % GM/WM/CSF files
        classFiles = {};
        for m = 1:3
            if m == 1, class = 'CSF'; elseif m == 2, class = 'GM'; elseif m == 3; class = 'WM'; end
            mask_regexp = ['^' b.subjs{nsub} '.*MNI.*' class '_probseg.nii.gz'];
            maskFile   = cellstr(spm_select('FPList', [b.derivDir b.subjs{nsub} '/anat/'], mask_regexp));
            if isempty(maskFile{1})
                error(['no file found for ',class,' !']);
            end
            classFiles{m} = maskFile;
            
            % unzip GM/WM/CSF - for now, conn requires .nii
            gunzip(maskFile);
        end
        
        % add to batch
        batch.Setup.masks.CSF.files{nsub}   = classFiles{1};
        batch.Setup.masks.Grey.files{nsub}  = classFiles{2};
        batch.Setup.masks.White.files{nsub} = classFiles{3};
        % specify that we only want 1 dimension
        batch.Setup.masks.CSF.dimensions   = 1;
        batch.Setup.masks.Grey.dimensions  = 1;
        batch.Setup.masks.White.dimensions = 1;
    end
    
    
    %% 4. GET CONFOUND REGRESSORS AND CONDITIONS
    
    fprintf('\nAdding regressors...\n\n');
    
    batch.Setup.conditions.names{1} = b.task; %only movie watching
    covar_names = {'motion','aCompCor','spikes'}; %from fmriprep_confound_regressors.m
    
    % add onsets and durations per session for target events, as well as
    % all covariates
    for nses = 1:b.nsessions
        for nsub = 1:length(b.subjs)
            % add task as whole session (would do same for resting state)
            batch.Setup.conditions.onsets{1}{nsub}{nses} = 0;
            batch.Setup.conditions.durations{1}{nsub}{nses} = inf;
            
            for c = 1:length(covar_names)
                c_name = covar_names{c};
                % grab covariates associated with this session for this subject:
                sessFiles = cellstr(spm_select('FPList', [b.covDir b.subjs{nsub} '/'], ['.*' b.task '_' c_name '.txt']));
                %add to batch
                if ~isempty(sessFiles{1}) %not all subjects will have a file for spike regressors
                    batch.Setup.covariates.names{c} = c_name;
                    batch.Setup.covariates.files{c}{nsub}{nses} = sessFiles{1};
                elseif isempty(sessFiles{1}) && ~strcmp(c_name,'spikes')
                    error(['no file found for ',c_name,' !']);
                end
            end
        end
    end
    
    % add confound covariates to denoising pipelines
    batch.Denoising.confounds.names = batch.Setup.covariates.names; %all (and only) my listed nuisance covariates
    for c = 1:length(batch.Denoising.confounds.names)
        batch.Denoising.confounds.deriv{c}=0; %do not add derviatives (note my motion files already have temporal derivs)
    end
    
    
    %% 5. DEFINE STEPS AND RUN
    
    % CONN Denoising --> detrending and regress out specified regressors from functional data
    batch.Denoising.filter     = [0.008,0.1];   %band-pass filter (Hz) - high pass and low pass
    batch.Denoising.detrending = 1;             %1: linear detrending, although redundant
    batch.Denoising.regbp      = 1;             %1= reg then bp filter - best if you have spikes. 2= Simultaneous regression and BP filtering
    
    % Additional setup steps
    batch.Setup.analyses        = 3; %just vox-to-vox, to visualize denoising, see conn_batch for codes
    batch.Setup.voxelmask       = 1; %1.Explicit mask -- note by default this would be SPM MNI brainmask.nii
    % here, I'm specifying the MNI 6 NLin Asym from templateflow (the template my data
    % was normalized to in fmriprep, in same resolution)
    batch.Setup.voxelmaskfile   = [base_dir 'templates/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz'];
    batch.Setup.voxelresolution = 3; %same as functional volumes (note that this is a code for conn, not the voxel size)
    batch.Setup.analysisunits   = 2; %2 = raw, 1 = PSC
    batch.Setup.outputfiles     = [0,1]; %(2) = create nifti for confound-corrected time series, see conn_batch
    batch.Setup.rois.names      = {'BrainMask'}; %note that this is just a filler given that we aren't running analyses
    batch.Setup.rois.files      = {batch.Setup.voxelmaskfile};
    
    
    batch.Denoising.done      = 1;
    batch.Denoising.overwrite = 1;
    batch.Analysis.done       = 0; % no analyses
    batch.vvAnalysis.done     = 0; % no analyses
    batch.dynAnalysis.done    = 0; % no analyses
    
    
    %%%%%%%%%%%%%%%%%%
    conn_batch(batch); % RUN!
    %%%%%%%%%%%%%%%%%%
    
    
    %% 6. Now move and rename the denoised .nii files for ease of use
    % Also converts zeros outside of brain mask to NaNs
    
    for nsub = 1:batch.Setup.nsubjects
        dataFile = fullfile(b.analysisDir,[b.subjs{nsub} '_task-' b.task],'results','preprocessing',['niftiDATA_Subject' num2str(nsub,'%03.f') '_Condition000.nii']);
        newFile  = fullfile(b.analysisDir,[b.subjs{nsub} '_task-' b.task '_MNI_denoised.nii']);
        
        % move file:
        unix(['mv ' dataFile ' ' newFile]);
        
        %convert file, zeros to NaNs (using cosmo - convenient for 4D data):
        data_new = cosmo_fmri_dataset(newFile);
        data_new.samples(data_new.samples == 0) = NaN;
        cosmo_map2fmri(data_new, newFile); %overwrite
    end
    
    
    %% 7. Remove unzipped niftis from fmriprep folders
    
    % fmriprep saves .nii.gz for space.
    % conn can accept .gz files, but it will automatically unzip them to work
    % with. Remove unzipped niftis in the fmriprep sub-directories:
    
    for nsub = 1:batch.Setup.nsubjects
        unix(['rm -rf ' b.derivDir b.subjs{nsub} '/*/*.nii']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end