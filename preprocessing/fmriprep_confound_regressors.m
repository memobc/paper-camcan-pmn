%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates confound regressor files for the CamCan movie watching dataset
% Jan 2020, Rose Cooper
%
% Loads tsv files generated from fmriprep to generate nuisance
% regressors per scan task/run/subject for data denoising
%
% saves 3 text files, with one column per confound:
% motion    --> FD, 6 realignment params and temporal derivatives (13 total)
% aCompCor  --> 6 PCs of a combined WM and CSF mask
% spikes    --> one regressors per outlying TR to effectively censor.
%               Criteria = FD > 0.6mm and.or STD DVARS > 2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;

b.scriptdir = pwd;
addpath(b.scriptdir);

addpath('/gsfs0/data/cooperrn/Documents'); %for tsvread()
base_dir = '/gsfs0/data/ritcheym/data/fmri/CamCan/pm-network/data/';

%where is my regressor information? (.tsv files from fmriprep)
b.derivDir  = [base_dir 'derivs/fmriprep/'];

task = 'movie';


%where to save new regressor files?
b.saveDir   = [base_dir 'derivs/confounds/'];
if ~exist(b.saveDir,'dir')
    mkdir(b.saveDir)
end

% load in file with *valid* subject IDs (subjects to include in analyses):
subjects = table2cell(readtable([base_dir 'valid_subjects_camcan-pmnetwork.csv'], 'ReadVariableNames', 0));



%% define confounds of interest

target_motion = {'framewise_displacement',...
    'trans_x','trans_x_derivative1',...
    'trans_y','trans_y_derivative1',...
    'trans_z','trans_z_derivative1',...
    'rot_x','rot_x_derivative1',...
    'rot_y','rot_y_derivative1',...
    'rot_z','rot_z_derivative1',...
    };
target_aComp = {'a_comp_cor_00','a_comp_cor_01','a_comp_cor_02',...
    'a_comp_cor_03','a_comp_cor_04','a_comp_cor_05',...
    };

spike_fd    = 0.6; % FD mm.
spike_dvars = 2;   % standardized variance in signal over voxels

fprintf('\nCreating confound regressors...\n');



%% loop through subjects

for i = 1:length(subjects)
    
    b.curSubj = subjects{i};
    %subject folder:
    b.subjDir   = [b.saveDir b.curSubj '/'];
    
    if ~exist(b.subjDir,'dir')
        mkdir(b.subjDir);
        fprintf('\n\t%s',b.curSubj);

        % fetch confound data:
        motionFile = [b.derivDir b.curSubj '/func/' b.curSubj '_task-' task '_desc-confounds_regressors.tsv'];
        [~,~,motionData] = tsvread(motionFile);
        %replace 'n\a' (first FD or first of derviatives) with 0
        motionData(2,contains(motionData(2,:),'n/a')) = {'0'};
        nTR = (size(motionData,1)-1);
        
        fprintf('\t%d timepoints',nTR);
        if nTR ~= 193
            error('Check number of scans');
        end
        
        
        %find target regressor columns for motion: ----------
        motion_Cols = cellfun(@(x) find(endsWith(motionData(1,:),x)),target_motion);
        %get target confound data:
        R = motionData(2:end,motion_Cols);
        %convert to double format
        R = cellstr(R);
        R = cellfun(@str2num, R);
        
        fileName = [b.subjDir b.curSubj '_' task '_motion.txt'];
        writetable(array2table(R),fileName,'Delimiter',' ','WriteVariableNames',false);
        
        
        %find target regressor columns for aCompCor: --------
        aComp_Cols = cellfun(@(x) find(endsWith(motionData(1,:),x)),target_aComp);
        %get target confound data:
        A = motionData(2:end,aComp_Cols);
        %convert to double format
        A = cellstr(A);
        A = cellfun(@str2num, A);
        
        fileName = [b.subjDir b.curSubj '_' task '_aCompCor.txt'];
        writetable(array2table(A),fileName,'Delimiter',' ','WriteVariableNames',false);
        
        
        % now create spike regressors --------------------
        spike_idx = [];
        
        % check to see if we need to add non-steady-state outliers (camcan
        % shouldn't have any, but just in case fmriprep thinks there are)
        nss = find(contains(motionData(1,:),'non_steady_state'));
        if ~isempty(nss)
            for n = 1:length(nss)
                idx = find(contains(motionData(2:end,nss(n)),'1.0'));
                spike_idx = [spike_idx idx];
            end
        end
        
        % now FD and dvars-based
        fd_col    = find(contains(motionData(1,:),'framewise_displacement'));
        dvars_col = find(contains(motionData(1,:),'std_dvars'));
        for tr = 1:size(R,1)
            if str2num(motionData{tr+1,fd_col}) > spike_fd || str2num(motionData{tr+1,dvars_col}) > spike_dvars
                spike_idx = [spike_idx tr];
            end
        end
        
        % save spikes file if there are any to exclude:
        if ~isempty(spike_idx)
            spike_regs = zeros(size(R,1),length(spike_idx));
            for s = 1:length(spike_idx)
                spike_regs(spike_idx(s),s) = 1;
            end
            fprintf('\t%d spikes found',size(spike_regs,2));
            
            fileName = [b.subjDir b.curSubj '_' task '_spikes.txt'];
            writetable(array2table(spike_regs),fileName,'Delimiter',' ','WriteVariableNames',false);
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end %end of loop through subjects ---------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%