%% Setup and Feature Extraction
% Ensure c1_dataVis.m has been run so filtered signals and labels exist.

% Parameters from Part 2.2b (Feature Selection)
WSize_sec = 0.1; 
Olap_pct = 0; 

% Group data for easy looping
signals_list = {filteredFlex, filteredPinch, filteredVF};
labels_list = {flexLabels, pinchLabels, VFLabels};
names = {'Flex', 'Pinch', 'VF'};

% Initialize struct to store the extracted features for each class
results = struct(); 

fprintf('Extracting features (WSize=%.2fs, Olap=%.2f)...\n', WSize_sec, Olap_pct);

for k = 1:3
    sig = signals_list{k};
    lbl = labels_list{k};
    
    % --- Windowing Logic ---
    WSize_samp = floor(WSize_sec * fs);
    nOlap = floor(Olap_pct * WSize_samp);
    hop = WSize_samp - nOlap;
    nx = length(sig);
    len = fix((nx - (WSize_samp - hop)) / hop);
    
    % Preallocate
    MAV_vec = zeros(1, len);
    VAR_vec = zeros(1, len);
    feat_lbl = zeros(1, len);
    
    % Get triggers for labeling
    Rise = gettrigger(lbl, 0.5);
    Fall = gettrigger(-lbl, -0.5);
    
    for i = 1:len
        idx_start = (i-1)*hop + 1;
        idx_end = idx_start + WSize_samp - 1;
        
        segment = sig(idx_start:idx_end);
        
        % Calculate Features
        MAV_vec(i) = mean(abs(segment));
        VAR_vec(i) = var(segment);
        
        % Labeling: 1 if window is strictly inside stimulation, 0 otherwise
        is_stim = any(idx_start >= Rise & idx_end <= Fall);
        feat_lbl(i) = double(is_stim);
    end
    
    % --- Separate Rest vs Stimulus Data ---
    % Store the features in the results struct
    results(k).Name = names{k};
    results(k).MAV_Stim = MAV_vec(feat_lbl == 1);
    results(k).MAV_Rest = MAV_vec(feat_lbl == 0);
    results(k).VAR_Stim = VAR_vec(feat_lbl == 1);
    results(k).VAR_Rest = VAR_vec(feat_lbl == 0);
end

%% Calculate Fisher's Discriminant Ratio (FDR)
% Formula: J = (mu1 - mu2)^2 / (var1^2 + var2^2)
% Note: Using Variance (sigma^2) directly in denominator

fprintf('\n=== Fisher''s Discriminant Ratio (FDR) Results ===\n');
fprintf('Higher value = Better discriminability\n');
fprintf('----------------------------------------------------------\n');
fprintf('%-25s | %-12s | %-12s\n', 'Comparison', 'FDR (MAV)', 'FDR (VAR)');
fprintf('----------------------------------------------------------\n');

% 1. Rest vs Stimulus (Internal comparison for each file)
for k = 1:3
    % Data for Stimulus vs Rest
    data_stim = results(k); 
    
    % --- MAV ---
    mu1 = mean(data_stim.MAV_Stim); var1 = var(data_stim.MAV_Stim);
    mu2 = mean(data_stim.MAV_Rest); var2 = var(data_stim.MAV_Rest);
    fdr_mav = ((mu1 - mu2)^2) / (var1 + var2);
    
    % --- VAR ---
    mu1 = mean(data_stim.VAR_Stim); var1 = var(data_stim.VAR_Stim);
    mu2 = mean(data_stim.VAR_Rest); var2 = var(data_stim.VAR_Rest);
    fdr_var = ((mu1 - mu2)^2) / (var1 + var2);
    
    fprintf('%-25s | %-12.4f | %-12.4f\n', [names{k} ' vs Rest'], fdr_mav, fdr_var);
end

fprintf('----------------------------------------------------------\n');

% 2. Stimulus vs Stimulus (Compare 'Stim' vector of one to 'Stim' vector of another)
pairs = [1 2; 1 3; 2 3]; % Flex-Pinch, Flex-VF, Pinch-VF

for p = 1:3
    idx1 = pairs(p, 1);
    idx2 = pairs(p, 2);
    
    name1 = names{idx1};
    name2 = names{idx2};
    
    % --- MAV ---
    mu1 = mean(results(idx1).MAV_Stim); var1 = var(results(idx1).MAV_Stim);
    mu2 = mean(results(idx2).MAV_Stim); var2 = var(results(idx2).MAV_Stim);
    fdr_mav = ((mu1 - mu2)^2) / (var1 + var2);
    
    % --- VAR ---
    mu1 = mean(results(idx1).VAR_Stim); var1 = var(results(idx1).VAR_Stim);
    mu2 = mean(results(idx2).VAR_Stim); var2 = var(results(idx2).VAR_Stim);
    fdr_var = ((mu1 - mu2)^2) / (var1 + var2);
    
    fprintf('%-25s | %-12.4f | %-12.4f\n', [name1 ' vs ' name2], fdr_mav, fdr_var);
end
fprintf('----------------------------------------------------------\n');