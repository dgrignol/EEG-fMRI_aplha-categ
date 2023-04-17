%% Toolboxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CosmoMVPA toolbox
run('/rds/projects/2018/hickeycm-insense/MATLAB_toolboxes/CoSMoMVPA-master/mvpa/cosmo_set_path.m')

% SPM toolbox (if needed)
% run('/rds/projects/2018/hickeycm-insense/EEG-fMRI/analys_scripts/MRI/initialise_spm.m')

% EEGlab toolbox
run('/rds/projects/2018/hickeycm-insense/MATLAB_toolboxes/eeglab14_1_2b/eeglab.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; 
close all;

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

    eeg_filepath = fullfile(root_dir,'Preproc_DATA/EEG/D0_Nt_Pd_200_250/',subj_str);

    % All clean trials, left hemi electrodes: ipsi; right: contra

    eeg_file_Dlat = [subj_str, '_EEG_Dlat.set']; % lateral distractor
    eeg_file_Tlat = [subj_str, '_EEG_Tlat.set']; % lateral target
    EEG_Dlat = pop_loadset('filename',eeg_file_Dlat,'filepath',eeg_filepath);
    EEG_Tlat = pop_loadset('filename',eeg_file_Tlat,'filepath',eeg_filepath);    
    EEG_Dlat = eeg_checkset(EEG_Dlat);
    EEG_Tlat = eeg_checkset(EEG_Tlat);
    
    if isubj==1
            nchan = EEG_Dlat.nbchan;
            Tlat_tfdata = {};
            Dlat_tfdata = {};
    end


    times = {};
    freqs = {};
    times2 = {};
    freqs2 = {};

    for electrode = 1:nchan 

        [ersp{electrode}, itc{electrode}, ...
            powbase{electrode}, times{electrode}, ...
            freqs{electrode}, erspboot{electrode}, ...
            itcboot{electrode}, tfdata_T{electrode}] = ... % , tfdata{1}{electrode}
            newtimef(EEG_Tlat.data(electrode,:,:), ...
            EEG_Tlat.pnts, [EEG_Tlat.xmin EEG_Tlat.xmax]*1000, EEG_Tlat.srate, [0.5 3], ...
            'plotphase', 'off', 'nfreqs', 26, 'timesout', 50, 'plotersp', 'off', 'plotitc', 'off', 'maxfreq', 52, 'baseline', 0,'trialbase', 'on', 'outputformat', 'old', 'scale', 'abs','verbose', 'off'); %'padratio', 16, 
        
        [ersp{electrode}, itc{electrode}, ...
            powbase{electrode}, times2{electrode}, ...
            freqs2{electrode}, erspboot{electrode}, ...
            itcboot{electrode}, tfdata_D{electrode}] = ... % , tfdata{1}{electrode}
            newtimef(EEG_Dlat.data(electrode,:,:), ...
            EEG_Dlat.pnts, [EEG_Dlat.xmin EEG_Dlat.xmax]*1000, EEG_Dlat.srate, [0.5 3], ...
            'plotphase', 'off', 'nfreqs', 26, 'timesout', 50, 'plotersp', 'off', 'plotitc', 'off', 'maxfreq', 52, 'baseline', 0,'trialbase', 'on', 'outputformat', 'old', 'scale', 'abs','verbose', 'off'); %'padratio', 16, 

        Tlat_tfdata{isubj, electrode} = tfdata_T{electrode};
        Dlat_tfdata{isubj, electrode} = tfdata_D{electrode};

    end
    Tlat_tfdata_subj = Tlat_tfdata{isubj};
    Dlat_tfdata_subj = Dlat_tfdata{isubj};

    save([save_dir '/P' sprintf('%01d',isubj) '_tf.mat'], ...
        'times', 'freqs',...
        'EEG_Tlat', 'EEG_Dlat', 'Tlat_tfdata_subj', 'Dlat_tfdata_subj','-v7.3')


end

save([save_dir '/' 'tf_ga.mat'], ...
      'times', 'freqs' ,'Tlat_tfdata','Dlat_tfdata','-v7.3')




