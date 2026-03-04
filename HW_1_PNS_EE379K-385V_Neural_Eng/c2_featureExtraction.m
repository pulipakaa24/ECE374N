%% Parameter Sweep Setup
% Note: Ensure you have run c1_dataVis.m first to get filteredFlex/fs/flexLabels
filteredSignal = filteredPinch; % Using Flex signal (High SNR) for demonstration 
label = pinchLabels;            % Labels of stimulus locations

% Sweep parameters
WSize_values = [0.05, 0.1, 0.3]; % Window sizes in seconds
Olap_values = [0, 0.25, 0.75];   % Overlap percentages

% Get trigger points once (indices in the raw signal)
Rise1 = gettrigger(label, 0.5);   % Start of stimulation
Fall1 = gettrigger(-label, -0.5); % End of stimulation

% Create Figure
figure('Name', 'Feature Extraction Sweep with MAV & VAR SNR', 'units', 'normalized', 'Position', [0, 0, 1, 1]);

plot_idx = 1; % Counter for subplot index

%% Loop through all combinations
for w = 1:length(WSize_values)
    for o = 1:length(Olap_values)
        
        % Current parameters
        current_WSize_s = WSize_values(w);
        current_Olap_pct = Olap_values(o);
        
        % Windowing calculations
        WSize_samp = floor(current_WSize_s * fs);      % Window size in samples
        nOlap = floor(current_Olap_pct * WSize_samp);  % Overlap in samples
        hop = WSize_samp - nOlap;                      % Hop size
        nx = length(filteredSignal);
        len = fix((nx - (WSize_samp - hop)) / hop);    % Total number of frames
        
        % Preallocate
        MAV_feature = zeros(1, len);
        VAR_feature = zeros(1, len);
        featureLabels = zeros(1, len);
        
        % Feature Extraction Loop
        for i = 1:len
            % Extract segment
            idx_start = (i-1)*hop + 1;
            idx_end = idx_start + WSize_samp - 1;
            
            % Check bounds
            if idx_end > nx
                break; 
            end
            
            segment = filteredSignal(idx_start:idx_end);
            
            % Calculate Features
            MAV_feature(i) = mean(abs(segment));
            VAR_feature(i) = var(segment);
            
            % Re-build label vector
            % Strict: Window must be fully inside stimulation to count as '1'
            is_stim = any(idx_start >= Rise1 & idx_end <= Fall1);
            featureLabels(i) = double(is_stim);
        end
        
        %% Calculate SNR
        % Separate Stimulus and Rest values using logical indexing
        stim_indices = (featureLabels == 1);
        rest_indices = (featureLabels == 0);
        
        % --- MAV SNR ---
        mean_mav_stim = mean(MAV_feature(stim_indices));
        mean_mav_rest = mean(MAV_feature(rest_indices));
        if mean_mav_rest > 0
            mav_snr = 20 * log10(mean_mav_stim / mean_mav_rest);
        else
            mav_snr = 0; 
        end
        
        % --- VAR SNR ---
        mean_var_stim = mean(VAR_feature(stim_indices));
        mean_var_rest = mean(VAR_feature(rest_indices));
        if mean_var_rest > 0
            var_snr = 20 * log10(mean_var_stim / mean_var_rest);
        else
            var_snr = 0; 
        end
        
        %% Plotting MAV (Left Column)
        subplot(9, 2, plot_idx);
        plot(MAV_feature, 'b', 'LineWidth', 0.5); hold on;
        plot(featureLabels * max(MAV_feature), 'r', 'LineWidth', 1);
        
        % Formatting Title with SNR
        title_str = sprintf('MAV (W=%.2g, Olap=%.2g) | SNR=%.2f dB', ...
                            current_WSize_s, current_Olap_pct, mav_snr);
        grid on;
        title(title_str, 'FontSize', 8); 
        
        if mod(plot_idx, 2) == 1
            ylabel('MAV');
        end
        xlim([1, len]);
        set(gca, 'XTickLabel', []); 
        
        plot_idx = plot_idx + 1;
        
        %% Plotting VAR (Right Column)
        subplot(9, 2, plot_idx);
        plot(VAR_feature, 'k', 'LineWidth', 0.5); hold on;
        plot(featureLabels * max(VAR_feature), 'r', 'LineWidth', 1);
        
        % Formatting VAR with SNR
        title_str_var = sprintf('VAR (W=%.2g, Olap=%.2g) | SNR=%.2f dB', ...
                                current_WSize_s, current_Olap_pct, var_snr);
        grid on;
        title(title_str_var, 'FontSize', 8);
        
        if mod(plot_idx, 2) == 0
            ylabel('VAR');
        end
        xlim([1, len]);
        set(gca, 'XTickLabel', []); 
        
        plot_idx = plot_idx + 1;
    end
end