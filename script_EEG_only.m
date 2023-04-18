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


    %% new EEG with swapped elec
    % new EEG swapped for left/right so that the diff btween new EEG
    % and orig EEG (ersp{3}) is contra vs ipsi on the left and 
    % ipsi vs contra on the right - no need of a baseline
    left = [EEG_Tlat.chanlocs(:).Y] > 0.0001;
    right = [EEG_Tlat.chanlocs(:).Y] < -0.0001;
    center = [EEG_Tlat.chanlocs(:).Y] > -0.0001 & [EEG_Tlat.chanlocs(:).Y] < 0.0001;
    sum(left) + sum(right) + sum(center)
    
    % swap Target lat eeg
    swap_EEG_Tlat = EEG_Tlat;
    swap_EEG_Tlat.data(left,:,:) = EEG_Tlat.data(right,:,:);
    swap_EEG_Tlat.data(right,:,:) = EEG_Tlat.data(left,:,:);
    % swap Distractor lat eeg
    swap_EEG_Dlat = EEG_Dlat;
    swap_EEG_Dlat.data(left,:,:) = EEG_Dlat.data(right,:,:);
    swap_EEG_Dlat.data(right,:,:) = EEG_Dlat.data(left,:,:);


    %% calculate ft
    if isubj==1
            nchan = EEG_Dlat.nbchan;
            cdn = zeros(nsubj, nchan,26,50);
            Dlat_ersp = zeros(nsubj, nchan,26,50);
            Tlat_ersp = zeros(nsubj, nchan,26,50);
            Tlat_tfdata = {};
            Dlat_tfdata = {};
            Tlat_tfdata_all = {};
    end

    ersp = {zeros(nchan,1)};
    itc = {zeros(nchan,1)};
    powbase = {zeros(nchan,1)};
    times = {zeros(nchan,1)};
    freqs = {zeros(nchan,1)};
    erspboot = {zeros(nchan,1)};
    itcboot = {zeros(nchan,1)};

    for electrode = 1:nchan 

        [ersp{electrode}, itc{electrode}, ...
            powbase{electrode}, times{electrode}, ...
            freqs{electrode}, erspboot{electrode}, ...
            itcboot{electrode}, tfdata{1}{electrode}] = ... % , tfdata{1}{electrode}
            newtimef({swap_EEG_Tlat.data(electrode,:,:) EEG_Tlat.data(electrode,:,:)}, ...
            EEG_Tlat.pnts, [EEG_Tlat.xmin EEG_Tlat.xmax]*1000, EEG_Tlat.srate, [0.5 3], ...
            'plotphase', 'off', 'timesout', 50, 'plotersp', 'off', 'plotitc', 'off', 'maxfreq', 52,'padratio', 8, 'baseline', NaN, 'outputformat', 'old', 'scale', 'abs','verbose', 'off'); %'padratio', 16, 

        cdn(isubj, electrode,:,:) = ersp{electrode}{3}; 
        Tlat_ersp(isubj, electrode,:,:) = ersp{electrode}{2};
        swap_Tlat_ersp(isubj, electrode,:,:) = ersp{electrode}{1};
        
%         Tlat_tfdata{electrode} = tfdata{electrode}{1};
%         Dlat_tfdata{electrode} = tfdata{electrode}{2};

%        Tlat_tfdata_all{isubj, electrode} = tfdata{electrode};
    end

%     save([save_dir '/P' sprintf('%01d',isubj) '_tf.mat'], ...
%         'ersp', 'itc', 'powbase', 'times', 'freqs', 'erspboot', 'itcboot', ...
%         'EEG_Tlat', 'EEG_Dlat', 'cdn', '-v7.3') %'Tlat_tfdata', 'Dlat_tfdata',


end

save([save_dir '/' 'tf_ga.mat'], ...
     'cdn', 'swap_Tlat_ersp', 'Tlat_ersp', 'times', 'freqs' ,'-v7.3') %'Tlat_tfdata_all',




