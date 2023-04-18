clear all; close all; 

eeglab nogui

working_directory = '\\its-rds.bham.ac.uk\rdsprojects\2018\hickeycm-insense\EEG-fMRI\alpha_power_analys\analys_scripts\';

output_directory = '\\its-rds.bham.ac.uk\rdsprojects\2018\hickeycm-insense\EEG-fMRI\alpha_power_analys\output_data\';



%% let's normalize:
if ~isfile([output_directory 'norm_tf_ga.mat'])
    
    % File does not exists. Create file.
    load([output_directory 'tf_26x50_baselineNaN_ersp_swapTlatVSTlat\' 'tf_ga.mat'])
    times = times{1};
    freqs = freqs{1}; 
    % load('example_tf.mat')
    % Tlat_tfdata = example;
    % clear Tlat_tfdata_all
    ciao = {Tlat_tfdata Dlat_tfdata};
%     % 1. absolute value
%     Tlat_tfdata = cellfun(@(x) abs(x).^2,Tlat_tfdata,'UniformOutput',false);
%     Dlat_tfdata = cellfun(@(x) abs(x).^2,Dlat_tfdata,'UniformOutput',false);
% 
%     % mean by trials
%     Tlat_tfdata = cellfun(@(x) squeeze(mean(x,3)),Tlat_tfdata,'UniformOutput',false);
%     Dlat_tfdata = cellfun(@(x) squeeze(mean(x,3)),Dlat_tfdata,'UniformOutput',false);
    
    for subj=1:25
        for chan=1:63
            %
            % 1. absolute value
            Tlat_tfdata{subj, chan} = abs(Tlat_tfdata{subj, chan}).^2; % should I divide by something?
            Dlat_tfdata{subj, chan} = abs(Dlat_tfdata{subj, chan}).^2; % {2} is Distractor

            % 2. normalize
            % traget
            size(Tlat_tfdata{subj, chan}) % freq times trials
            idx = find(times < 0, 1, 'last' ); % time 0
            P = Tlat_tfdata{subj, chan};
            Pbaseline = mean(P(:,1:idx,:),2); % time before baseline
            Pbaseline_mat = repmat(Pbaseline, [1, size(P,2),1]);
            Tlat_tfdata{subj, chan} = (P - Pbaseline_mat) ./ (Pbaseline_mat);
            size(Tlat_tfdata{subj, chan})
            % Distractor
            size(Dlat_tfdata{subj, chan}) % freq times trials
            idx = find(times < 0, 1, 'last' ); % time 0
            P = Dlat_tfdata{subj, chan};
            Pbaseline = mean(P(:,1:idx,:),2); % time before baseline
            Pbaseline_mat = repmat(Pbaseline, [1, size(P,2),1]);
            Dlat_tfdata{subj, chan} = (P - Pbaseline_mat) ./ (Pbaseline_mat);
            size(Dlat_tfdata{subj, chan})

            % mean by trials
            Tlat_tfdata{subj, chan} = squeeze(mean(Tlat_tfdata{subj, chan},3));
            Dlat_tfdata{subj, chan} = squeeze(mean(Dlat_tfdata{subj, chan},3));

            % put in a matrix with following dimensions: subj/chan/freq/time
            Tlat_mat(subj,chan,:,:) = Tlat_tfdata{subj, chan};
            Dlat_mat(subj,chan,:,:) = Dlat_tfdata{subj, chan}; 
            
        end
    end
    save([output_directory 'norm_tf_ga.mat'], ...
        'Tlat_mat', 'Dlat_mat', 'times', 'freqs', ...
        '-v7.3')
end
% File exist. Load file.
load([output_directory 'norm_tf_ga.mat'])
     


%% select electrodes
% left hemi electrodes: ipsi; right: contra
% p5/6 po7/8 p7/8
load('chanlocs.mat')

load('chanlabels')
chanlabels = {chanlabels};
list_of_chan_ipsi =  {'P5' 'PO7' 'P7'}; 
idx_ipsi = zeros(length(list_of_chan_ipsi),1);

list_of_chan_contra =  {'P6' 'PO8' 'P8'}; 
idx_contra = zeros(length(list_of_chan_contra),1);

for i=1:length(list_of_chan_contra)
    
    idx_ipsi(i) = find(ismember(chanlabels{:}, list_of_chan_ipsi{i}));
    idx_contra(i) = find(ismember(chanlabels{:}, list_of_chan_contra{i}));
end
% po7 = [chanlocs(:).Y] > 1;
% po8 = [chanlocs(:).Y] < -1;


%% Target plots
% reduce electrodes to 3 contra & 3 ipsi 
contra = Tlat_mat(:,idx_contra,:,:); % subj/chan/freq/time
ipsi = Tlat_mat(:,idx_ipsi,:,:);
size(ipsi)

% difference contra - ipsi
diff = contra - ipsi;
contra = squeeze(mean(contra,2));
ipsi = squeeze(mean(ipsi,2));
diff = squeeze(mean(diff,2));
size(contra)
size(diff)

% subj mean for contra, ipsi & difference
size(contra)
contra = squeeze(mean(contra,1)); % mean subj contra
size(contra)

size(ipsi)
ipsi = squeeze(mean(ipsi,1)); % mean subj ipsi
size(ipsi)

size(diff)
diff = squeeze(mean(diff,1)); % mean subj difference
size(diff)

% plot
figure; 
subplot(1,3,1)
tftopo(contra, times, freqs);
title('contra')
subplot(1,3,2)
tftopo(ipsi, times, freqs);
title('ipsi')
subplot(1,3,3)
tftopo(diff, times, freqs);
title('contra - ipsi')
suptitle(['Target Lateral' newline 'Chan: P5/6 PO7/8 P7/8'])

%topoplot
min_freq = 7;
max_freq = 12;
% % determine plots high and low limits
% a = Tlat_mat(:,:,freqs>=min_freq & freqs<=max_freq,:);
% a = squeeze(mean(a,3));
% a = squeeze(mean(a,3));
% a = squeeze(mean(a,1));
% lo_hi = [min(a) max(a)]; %topoplot value limits
lo_hi = [-0.5 1.3]
figure;
time_windows = [-128 0; 0 50; 50 100; 100 150; 150 200; 200 250; 250 300; 300 350; 350 428]; 
for i=1:9
    
    subplot(3,3,i)
    
    %     min_time = 428/8*(i-1);
    %     max_time = 428/8*(i);
    min_time = time_windows(i,1);
    max_time = time_windows(i,2);    
    
    title(['[' num2str(min_time) ' ' num2str(max_time) ']'])
    
    a = Tlat_mat(:,:,freqs>=min_freq & freqs<=max_freq,times>=min_time & times<=max_time);
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,1));

    % left hemi electrodes: ipsi; right: contra
    Tlat_topop = squeeze(mean(a,[1 3 4]));
    size(Tlat_topop)
    topoplot(Tlat_topop,chanlocs,'maplimits',lo_hi)

end
suptitle(['Target Lateral' newline 'left: ipsi - - - right: contra'])

%% Distractor plots
% reduce electrodes to 3 contra & 3 ipsi 
contra = Dlat_mat(:,idx_contra,:,:); % subj/chan/freq/time
ipsi = Dlat_mat(:,idx_ipsi,:,:);
size(ipsi)

% difference contra - ipsi
diff = contra - ipsi;
contra = squeeze(mean(contra,2));
ipsi = squeeze(mean(ipsi,2));
diff = squeeze(mean(diff,2));
size(contra)
size(diff)

% subj mean for contra, ipsi & difference
size(contra)
contra = squeeze(mean(contra,1)); % mean subj contra
size(contra)

size(ipsi)
ipsi = squeeze(mean(ipsi,1)); % mean subj ipsi
size(ipsi)

size(diff)
diff = squeeze(mean(diff,1)); % mean subj difference
size(diff)

% plot
figure; 
subplot(1,3,1)
tftopo(contra, times, freqs);
title('contra')
subplot(1,3,2)
tftopo(ipsi, times, freqs);
title('ipsi')
subplot(1,3,3)
tftopo(diff, times, freqs);
title('contra - ipsi')
suptitle(['Distractor Lateral' newline 'Chan: P5/6 PO7/8 P7/8'])

%topoplot
% all data Dlat_mat
figure;
min_freq = 8;
max_freq = 12;
time_windows = [-128 0; 0 50; 50 100; 100 150; 150 200; 200 250; 250 300; 300 350; 350 428]; 
for i=1:9
    
    subplot(3,3,i)
    
    %     min_time = 428/8*(i-1);
    %     max_time = 428/8*(i);
    min_time = time_windows(i,1);
    max_time = time_windows(i,2);    
    
    title(['[' num2str(min_time) ' ' num2str(max_time) ']'])
    
    a = Dlat_mat(:,:,freqs>7 & freqs<13,times>=min_time & times<=max_time);
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,1));

    % left hemi electrodes: ipsi; right: contra
    Dlat_topop = squeeze(mean(a,[1 3 4]));
    size(Dlat_topop)
%     'maplimits'       
%       'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
%       'maxmin'   -> scale colors to the data range (makes green mid-range); 
    topoplot(Dlat_topop,chanlocs,'maplimits','maxmin')

end
suptitle(['Distractor Lateral' newline 'left: ipsi - - - right: contra'])


%% Target vs Distractor plots
TDlat_mat = Tlat_mat - Dlat_mat;

% reduce electrodes to 3 contra & 3 ipsi 
contra = TDlat_mat(:,idx_contra,:,:); % subj/chan/freq/time
ipsi = TDlat_mat(:,idx_ipsi,:,:);
size(ipsi)

% difference contra - ipsi
diff = contra - ipsi;
contra = squeeze(mean(contra,2));
ipsi = squeeze(mean(ipsi,2));
diff = squeeze(mean(diff,2));
size(contra)
size(diff)

% subj mean for contra, ipsi & difference
size(contra)
contra = squeeze(mean(contra,1)); % mean subj contra
size(contra)

size(ipsi)
ipsi = squeeze(mean(ipsi,1)); % mean subj ipsi
size(ipsi)

size(diff)
diff = squeeze(mean(diff,1)); % mean subj difference
size(diff)

% plot
figure; 
subplot(1,3,1)
tftopo(contra, times, freqs);
title('contra')
subplot(1,3,2)
tftopo(ipsi, times, freqs);
title({'','ipsi'})
subplot(1,3,3)
tftopo(diff, times, freqs);
title('contra - ipsi')
suptitle(['Target - Distractor Lateral' newline 'Chan: P5/6 PO7/8 P7/8' ])

%topoplot
% all data TDlat_mat
figure;
min_freq = 8;
max_freq = 12;
time_windows = [-128 0; 0 50; 50 100; 100 150; 150 200; 200 250; 250 300; 300 350; 350 428]; 
for i=1:9
    
    subplot(3,3,i)
    
    %     min_time = 428/8*(i-1);
    %     max_time = 428/8*(i);
    min_time = time_windows(i,1);
    max_time = time_windows(i,2);    
    
    title(['[' num2str(min_time) ' ' num2str(max_time) ']'])
    
    a = TDlat_mat(:,:,freqs>7 & freqs<13,times>=min_time & times<=max_time);
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,3));
    a = squeeze(mean(a,1));

    % left hemi electrodes: ipsi; right: contra
    TDlat_topop = squeeze(mean(a,[1 3 4]));
    size(TDlat_topop)
%     'maplimits'       
%       'absmax'   -> scale map colors to +/- the absolute-max (makes green 0); 
%       'maxmin'   -> scale colors to the data range (makes green mid-range); 
    topoplot(TDlat_topop,chanlocs,'maplimits','maxmin')

end
suptitle(['Target - Distractor Lateral' newline 'left: ipsi - - - right: contra'])