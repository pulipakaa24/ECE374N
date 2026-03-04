%% 
clear all
close all
clc

%% Import Data
load('data.mat');

%% Example: Plot the raw signal
flexSignal=Flex.signal;
pinchSignal = Pinch.signal;
VFSignal = VF.signal;
flexLabels=Flex.trigger;
pinchLabels=Pinch.trigger;
VFLabels = VF.trigger;

%% Plot the raw signals
% Grouping variables into cell arrays to loop through them cleanly
signals = {flexSignal, pinchSignal, VFSignal};
labels_all = {flexLabels, pinchLabels, VFLabels};
titles = {'Raw Flex Signal with Stimulation Pattern (Yellow)', ...
          'Raw Pinch Signal with Stimulation Pattern (Yellow)', ...
          'Raw VF Signal with Stimulation Pattern (Yellow)'};

% Create one large figure for all three subplots
figure('units','normalized','Position',[0.1, 0.1, 0.7, 0.8])

for i = 1:3
    sig = signals{i};
    lbls = labels_all{i};
    
    % Find trigger start and end points using your custom function
    TRIG = gettrigger(lbls, 0.5);
    TRIGend = gettrigger(-lbls, -0.5);
    
    % Create a subplot (3 rows, 1 column, current index i)
    subplot(3, 1, i);
    
    % Plot normalized signal and labels
    plot((1:length(sig))./fs, zscore(sig));
    hold on;
    plot((1:length(sig))./fs, zscore(lbls), 'y', 'LineWidth', 1.5);
    
    % Plot stem markers for triggers (added conditional checks just in case a signal has no triggers)
    if ~isempty(TRIG)
        stem(TRIG./fs, ones(length(TRIG),1)*max(zscore(lbls)), 'Color', 'g');
    end
    if ~isempty(TRIGend)
        stem(TRIGend./fs, ones(length(TRIGend),1)*max(zscore(lbls)), 'Color', 'r');
    end
    
    % Formatting
    grid on; grid minor;
    xlim([0, length(sig)./fs]);
    xlabel('Time (s)');
    % Note: Because you used zscore(), the amplitude is no longer in uV, but in standard deviations
    ylabel('Amplitude (stdDevs)'); 
    title(titles{i});
end

%% Example: PSD estimates
figure('Name', 'Raw PSD Estimates', 'units','normalized','Position',[0.1,0.1,0.5,0.5])
h = spectrum.welch; 

% --- REST (All signals combined) ---
% Find indices where labels are 0 for each signal type
[flex_rows_rest, ~, ~]  = find(flexLabels == 0);
[pinch_rows_rest, ~, ~] = find(pinchLabels == 0);
[vf_rows_rest, ~, ~]    = find(VFLabels == 0);

% Extract the actual signal data during those rest periods
flex_restData  = flexSignal(flex_rows_rest);
pinch_restData = pinchSignal(pinch_rows_rest);
vf_restData    = VFSignal(vf_rows_rest);

% Concatenate all rest data into one large vector for a robust estimate
all_restData = [flex_restData; pinch_restData; vf_restData];

% Calculate and Plot Rest PSD in Black
SOIf_rest = psd(h, all_restData, 'Fs', fs); 
plot(SOIf_rest.Frequencies, 10*log10(SOIf_rest.Data), 'k', 'LineWidth', 1.5); 
hold on;

% --- FLEX ---
[flex_rows_act, ~, ~] = find(flexLabels>0);
flex_signalOfInterest = flexSignal(flex_rows_act);
SOIf_flex = psd(h, flex_signalOfInterest, 'Fs', fs); 
plot(SOIf_flex.Frequencies, 10*log10(SOIf_flex.Data), 'g'); % Plot in Green

% --- PINCH ---
[pinch_rows_act, ~, ~] = find(pinchLabels>0);
pinch_signalOfInterest = pinchSignal(pinch_rows_act);
SOIf_pinch = psd(h, pinch_signalOfInterest, 'Fs', fs); 
plot(SOIf_pinch.Frequencies, 10*log10(SOIf_pinch.Data), 'r'); % Plot in Red

% --- VF ---
[VF_rows_act, ~, ~] = find(VFLabels>0);
VF_signalOfInterest = VFSignal(VF_rows_act);
SOIf_VF = psd(h, VF_signalOfInterest, 'Fs', fs); 
plot(SOIf_VF.Frequencies, 10*log10(SOIf_VF.Data), 'b'); % Plot in Blue

% --- FORMATTING ---
grid on;
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Rest', 'Flex', 'Pinch', 'VF');
title('Power Spectral Density Estimates (Raw Signals)');

%% Bandpass Filtering
% 
fc1 = 800; % first cutoff frequency in Hz 
fc2 = 2200; % second cutoff frequency in Hz 
% 
% % normalize the frequencies
Wp = [fc1 fc2]*2/fs;

% Build a Butterworth bandpass filter of 4th order
% check the "butter" function in matlab
n = 2; 
[b, a] = butter(n, Wp, 'bandpass');

% Filter data of both classes with a non-causal filter
% Hint: use "filtfilt" function in MATLAB
% filteredSignal = ;

filteredFlex  = filtfilt(b, a, flexSignal);
filteredPinch = filtfilt(b, a, pinchSignal);
filteredVF    = filtfilt(b, a, VFSignal);

%% Compare VF Signal Before and After Filtering (Time Domain)

% Group the raw and filtered VF signals to loop through them cleanly
vf_signals = {VFSignal, filteredVF};
vf_labels = {VFLabels, VFLabels}; % Labels remain the exact same
vf_titles = {'Raw VF Signal with Stimulation Pattern (Yellow)', ...
             'Filtered VF Signal (800-2200 Hz) with Stimulation Pattern (Yellow)'};

% Create one figure for the two subplots (slightly shorter since it's only 2 plots)
figure('Name', 'VF Filter Comparison (Time Domain)', 'units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.6])

for i = 1:2
    sig = vf_signals{i};
    lbls = vf_labels{i};
    
    % Find trigger start and end points using your custom function
    TRIG = gettrigger(lbls, 0.5);
    TRIGend = gettrigger(-lbls, -0.5);
    
    % Create a subplot (2 rows, 1 column, current index i)
    subplot(2, 1, i);
    
    % Plot normalized signal and labels
    plot((1:length(sig))./fs, zscore(sig));
    hold on;
    plot((1:length(sig))./fs, zscore(lbls), 'y', 'LineWidth', 1.5);
    
    % Plot stem markers for triggers
    if ~isempty(TRIG)
        stem(TRIG./fs, ones(length(TRIG),1)*max(zscore(lbls)), 'Color', 'g');
    end
    if ~isempty(TRIGend)
        stem(TRIGend./fs, ones(length(TRIGend),1)*max(zscore(lbls)), 'Color', 'r');
    end
    
    % Formatting
    grid on; grid minor;
    xlim([0, length(sig)./fs]);
    xlabel('Time (s)');
    ylabel('Amplitude (stdDevs)'); 
    title(vf_titles{i});
end


%% Compare VF Signal Before and After Filtering (Frequency Domain / PSD)

figure('Name', 'VF PSD Comparison', 'units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.4])

% --- RAW VF ---
[VF_rows_act, ~, ~] = find(VFLabels>0);
VF_raw_signalOfInterest = VFSignal(VF_rows_act);
h = spectrum.welch;
SOIf_VF_raw = psd(h, VF_raw_signalOfInterest, 'Fs', fs); 
plot(SOIf_VF_raw.Frequencies, 10*log10(SOIf_VF_raw.Data), 'b'); % Plot Raw in Blue
hold on; 

% --- FILTERED VF ---
% Reuse the exact same trigger rows since the timing hasn't changed
VF_filt_signalOfInterest = filteredVF(VF_rows_act);
SOIf_VF_filt = psd(h, VF_filt_signalOfInterest, 'Fs', fs); 
plot(SOIf_VF_filt.Frequencies, 10*log10(SOIf_VF_filt.Data), 'r'); % Plot Filtered in Red

% --- FORMATTING ---
grid on; grid minor;
title('PSD of VF Signal Before and After Filtering');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Raw VF', 'Filtered VF (800 - 2200 Hz)');
% xlim([0, 3000]); % Zoom in to see the filter cutoff roll-off clearly