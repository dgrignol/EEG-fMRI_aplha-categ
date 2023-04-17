%% Toolboxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CosmoMVPA toolbox
run('/rds/projects/2018/hickeycm-insense/MATLAB_toolboxes/CoSMoMVPA-master/mvpa/cosmo_set_path.m')

% SPM toolbox (if needed)
% run('/rds/projects/2018/hickeycm-insense/EEG-fMRI/analys_scripts/MRI/initialise_spm.m')

% EEGlab toolbox
run('/rds/projects/2018/hickeycm-insense/MATLAB_toolboxes/eeglab14_1_2b/eeglab.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Participants numbers after cleaning
subj_num = [1:22,24,29:30];
nsubj = length(subj_num);

root_dir = '/rds/projects/2018/hickeycm-insense/EEG-fMRI';

% Save data in a subdirectory of:
save_dir = fullfile(root_dir,'alpha_power_analys','output_data');

% Common grey matter mask for all participants (searchlight or any area of the brain)
mask_file = fullfile(root_dir,'DATA','group_mask_grey.nii');

% If you want to work on the OSC only, you can use this mask 
% (generated from Target classification searchlight):
% mask_file = fullfile(root_dir,'results','MRI','MVPA2_2cat_GLMsingle_searchlights','results','H_accT.nii');

for isubj = 1:nsubj

    subj_str = ['P' sprintf('%02.f',subj_num(isubj))];


    %% Loading behavioural data (trial info for all trials)
    all_trials = [];
    for iblock = 1:4
        beh_tbl = readtable(fullfile(root_dir,'DATA', subj_str, 'behav',[subj_str '_block' num2str(iblock) '.csv']));
        all_trials = [all_trials; beh_tbl];
    end
    % From this table, you can extract any information about any trial.

    %% Loading epoched EEG data

    eeg_filepath = fullfile(root_dir,'Preproc_DATA/EEG/D0_Nt_Pd/',subj_str);

    % All clean trials, left hemi electrodes: ipsi; right: contra

    eeg_file_Dlat = [subj_str, '_EEG_Dlat.set']; % lateral distractor
    eeg_file_Tlat = [subj_str, '_EEG_Tlat.set']; % lateral target
    EEG_Dlat = pop_loadset('filename',eeg_file_Dlat,'filepath',eeg_filepath);
    EEG_Tlat = pop_loadset('filename',eeg_file_Tlat,'filepath',eeg_filepath);    
    EEG_Dlat = eeg_checkset(EEG_Dlat);
    EEG_Tlat = eeg_checkset(EEG_Tlat);

    % Block and trial number for each epoch are in the event at latency 0
    % ms. For example: B 1 T 3 (for block 1, trial 3).
    % To get the indices from 1 to 512 from the epoched EEG data, 
    % use the get_trial_idx_list() function:
    trials_Dlat = get_trial_idx_list(EEG_Dlat);
    trials_Tlat = get_trial_idx_list(EEG_Tlat);

    % Target category of T lateral EEG trials
    Tcat_Tlat = all_trials.catT(trials_Tlat);

    % If you want to use the eeglab GUI:
    % ALLEEG = [];
    % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG_Dlat, 0 );
    % [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG_Tlat, 0 );
    % eeglab redraw

    % From EEG_Tlat and EEG_Dlat, you can calculate trial-wise alpha power 
    % or any other measure.


    %% Loading corresponding MRI data in CosmoMVPA (beta weights for each trial)
    
    fmri_subj_folder = fullfile(root_dir,'Preproc_DATA/MRI/Attention/betas_from_GLMsingle', subj_str,'all_trials');
    
    % Say we want all trials with a lateral target:
    % (if we want all trials, just load trials from 1 to 512 instead of using trials_Tlat here)
    nTlat = length(trials_Tlat);
    ds_fmri = cell(nTlat,1);

    for itrial = 1:nTlat
        fmri_filename = [subj_str '_Trial' num2str(trials_Tlat(itrial)) '_beta.nii'];
        ds_fmri{itrial} = cosmo_fmri_dataset(fullfile(fmri_subj_folder,fmri_filename),'targets',1,'chunks',subj_num(isubj),'mask',mask_file);
    end
    ds_fmri = cosmo_stack(ds_fmri);

    % Use data from "all_trials" to populate ds_fmri.sa.targets and
    % ds_fmri.sa.chunks, depending on what you want to do.
    ds_fmri.sa.targets = Tcat_Tlat;
    cosmo_check_dataset(ds_fmri);
    % To remove trials: cosmo_slice
    % To average trials: cosmo_average_samples

    % To prepare train-test sets:
    % cosmo_chunkize
    % cosmo_nchoosek_partitioner
    

    % Then you run your classification

    % See /rds/projects/2018/hickeycm-insense/EEG-fMRI/analys_scripts/EEG_MRI_2
    % for examples of how to do it

    % Then you can run your EEG-fMRI correlation across trials

    % To save data from CosmoMVPA into .nii file, use:
    % cosmo_map2fmri()


end


% cosmo_stat
% cosmo_montecarlo_cluster_stat
