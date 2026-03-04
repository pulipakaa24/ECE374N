%% c3_classification_complete.m
% This script performs 10-fold cross-validation for Part 2.3 of the assignment.
% It answers:
% 1. Classifiability of Stimuli vs Rest
% 2. Classifiability of Stimuli vs Stimuli
% 3. Comparison of MAV vs VAR features
% 4. Evaluation of Confusion Matrices

clearvars -except filteredFlex filteredPinch filteredVF flexLabels pinchLabels VFLabels fs;
clc;

% Check if data is loaded
if ~exist('filteredFlex', 'var')
    error('Error: Filtered signals not found. Please run c1_dataVis.m first.');
end

%% 1. Feature Extraction (WSize=100ms, Olap=0)
fprintf('1. Extracting Features (100ms Window, 0%% Overlap)...\n');

WSize_sec = 0.1; 
Olap_pct = 0; 
WSize = floor(WSize_sec * fs);
nOlap = floor(Olap_pct * WSize);
hop = WSize - nOlap;

% Organize data for looping
% Index 1=VF, 2=Flex, 3=Pinch
sigs = {filteredVF, filteredFlex, filteredPinch};
lbls = {VFLabels, flexLabels, pinchLabels};
names = {'VF', 'Flex', 'Pinch'};

feats = struct(); % Structure to hold features

for k = 1:3
    sig = sigs{k};
    lab = lbls{k};
    
    nx = length(sig);
    len = fix((nx - (WSize - hop)) / hop);
    
    MAV_vec = zeros(1, len);
    VAR_vec = zeros(1, len);
    LBL_vec = zeros(1, len);
    
    Rise = gettrigger(lab, 0.5);
    Fall = gettrigger(-lab, -0.5);
    
    for i = 1:len
        idx_start = (i-1)*hop + 1;
        idx_end = idx_start + WSize - 1;
        segment = sig(idx_start:idx_end);
        
        MAV_vec(i) = mean(abs(segment));
        VAR_vec(i) = var(segment);
        
        % Label: 1 if window is strictly inside stimulation
        is_stim = any(idx_start >= Rise & idx_end <= Fall);
        LBL_vec(i) = double(is_stim);
    end
    
    feats(k).MAV = MAV_vec;
    feats(k).VAR = VAR_vec;
    feats(k).LBL = LBL_vec;
    feats(k).Name = names{k};
end

%% 2. Define Comparisons
% We need to run classification for these specific pairs:
comparisons = {
    'VF vs Rest',      1, 0;  % 0 denotes "Rest" class
    'Flex vs Rest',    2, 0;
    'Pinch vs Rest',   3, 0;
    'Flex vs Pinch',   2, 3;
    'Flex vs VF',      2, 1;
    'Pinch vs VF',     3, 1;
};

%% 3. Classification Loop (10-Fold CV)
fprintf('\n2. Running 10-Fold Cross Validation...\n');
fprintf('----------------------------------------------------------------\n');
fprintf('%-20s | %-12s | %-12s | %-15s\n', 'Comparison', 'Acc (MAV)', 'Acc (VAR)', 'Best Feature');
fprintf('----------------------------------------------------------------\n');

for c = 1:size(comparisons, 1)
    comp_name = comparisons{c, 1};
    idx1 = comparisons{c, 2};
    idx2 = comparisons{c, 3};
    
    % --- Prepare Data for Class 1 ---
    % Get Stimulus features (Label == 1)
    f1_MAV = feats(idx1).MAV(feats(idx1).LBL == 1);
    f1_VAR = feats(idx1).VAR(feats(idx1).LBL == 1);
    
    % --- Prepare Data for Class 2 (or Rest) ---
    if idx2 == 0 
        % If comparing vs Rest, get Rest features (Label == 0) from the SAME signal
        f2_MAV = feats(idx1).MAV(feats(idx1).LBL == 0);
        f2_VAR = feats(idx1).VAR(feats(idx1).LBL == 0);
        label_names = {feats(idx1).Name, 'Rest'};
    else
        % If comparing vs another Stimulus, get Stimulus features (Label == 1)
        f2_MAV = feats(idx2).MAV(feats(idx2).LBL == 1);
        f2_VAR = feats(idx2).VAR(feats(idx2).LBL == 1);
        label_names = {feats(idx1).Name, feats(idx2).Name};
    end
    
    % Combine Data
    X_MAV = [f1_MAV, f2_MAV]'; % Transpose to column vector
    X_VAR = [f1_VAR, f2_VAR]';
    
    % Create Labels (1 for Class 1, 2 for Class 2)
    Y = [ones(length(f1_MAV), 1); 2 * ones(length(f2_MAV), 1)];
    
    % --- 10-Fold Cross Validation ---
    k = 10;
    cv = cvpartition(Y, 'KFold', k); % Random split (answering Q4!)
    
    acc_mav = 0;
    acc_var = 0;
    conf_mav = zeros(2,2); % Accumulate confusion matrix
    
    for i = 1:k
        train_idx = cv.training(i);
        test_idx = cv.test(i);
        
        % MAV Classification
        pred_mav = classify(X_MAV(test_idx), X_MAV(train_idx), Y(train_idx));
        acc_mav = acc_mav + sum(pred_mav == Y(test_idx)) / length(pred_mav);
        
        % Build Confusion Matrix for MAV (just one example needed for assignment)
        % Rows = True Class, Cols = Predicted Class
        current_conf = confusionmat(Y(test_idx), pred_mav);
        % Handle edge case if a fold misses a class
        if size(current_conf,1) == 2
             conf_mav = conf_mav + current_conf;
        end

        % VAR Classification
        pred_var = classify(X_VAR(test_idx), X_VAR(train_idx), Y(train_idx));
        acc_var = acc_var + sum(pred_var == Y(test_idx)) / length(pred_var);
    end
    
    % Average Accuracy
    mean_acc_mav = (acc_mav / k) * 100;
    mean_acc_var = (acc_var / k) * 100;
    
    % Determine Winner
    if mean_acc_mav > mean_acc_var
        winner = 'MAV';
    elseif mean_acc_var > mean_acc_mav
        winner = 'VAR';
    else
        winner = 'Tie';
    end
    
    fprintf('%-20s | %-11.1f%% | %-11.1f%% | %-15s\n', comp_name, mean_acc_mav, mean_acc_var, winner);
    
    % --- Display Confusion Matrix for MAV ---
    % Only printing logic to keep output clean, answering "Observe confusion matrices"
    fprintf('   Confusion Matrix (MAV) for %s:\n', comp_name);
    fprintf('   True %-6s: [ %4d  %4d ] (Predicted %s / %s)\n', label_names{1}, conf_mav(1,1), conf_mav(1,2), label_names{1}, label_names{2});
    fprintf('   True %-6s: [ %4d  %4d ]\n\n', label_names{2}, conf_mav(2,1), conf_mav(2,2));
end
fprintf('----------------------------------------------------------------\n');

%% 4. Answer Prompts
fprintf('\n=== Automated Analysis for Part 2.3 ===\n');
fprintf('1. Check "Pinch vs Rest" accuracy above. Is it low? (Likely yes, due to low SNR).\n');
fprintf('2. Check "Flex vs Pinch". Can they be distinguished?\n');
fprintf('3. Observe the Confusion Matrices: Are they balanced? \n');
fprintf('   - If one class is predicted much more often, the classifier is biased.\n');
fprintf('4. Feature Performance: Look at the "Best Feature" column.\n');
fprintf('   - MAV is typically more robust for these signals.\n');
fprintf('5. Validation Fairness (Assignment Q4):\n');
fprintf('   - This script uses "cvpartition", which splits data RANDOMLY.\n');
fprintf('   - Since EMG/ENG signals are time-series, random splitting causes "data leakage"\n');
fprintf('     (training on samples immediately adjacent to test samples).\n');
fprintf('   - Therefore, this is likely NOT a fair assessment of generalization.\n');