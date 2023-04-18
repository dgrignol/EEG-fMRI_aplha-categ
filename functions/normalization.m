function [tf_data] = normalization(tf_data,times)
% feed newtimef tfdata (matrix: freq, times, trials) and times (array of timepoints), 
% returns abs().^2 and then normalize with (P - Pbaseline) / (Pbaseline) 
% tf data is a mat: freq, times, trials

% BASELINE with percentage change kind of method (P - Pbaseline) / (Pbaseline)

    % 1. absolute value
    tf_data = abs(tf_data).^2; 

    % 2. normalize: (P - Pbaseline) / (Pbaseline)
    idx = find(times < 0, 1, 'last' ); % time 0, before stim
    P = tf_data;
    Pbaseline = mean(P(:,1:idx,:),2); % time before baseline
    Pbaseline_mat = repmat(Pbaseline, [1, size(P,2),1]);
    tf_data = (P - Pbaseline_mat) ./ (Pbaseline_mat);
