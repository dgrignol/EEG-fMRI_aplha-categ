function [kept_trials] = get_trial_idx_list(EEG)

% This function takes as an input the EEGlab EEG structure, searches events
% with name 'BX_TY' and returns the index (from 1 to 512) of the trials
% that are in the structure

% Creating the full list of trigger types that index block and trial number
trial_idx = cell(512,1);
idx = 1;
for iblock = 1:4
for itrial = 1:128
    trial_idx{idx} = ['B' num2str(iblock) '_T' num2str(itrial)];
    idx=idx+1;
end
end

kept_trials = [];
for ievent = 1:length(EEG.event)
    if strcmp(EEG.event(ievent).type(1),'B')
        kept_trials = [kept_trials; find(strcmp(EEG.event(ievent).type,trial_idx))];
    end
end
