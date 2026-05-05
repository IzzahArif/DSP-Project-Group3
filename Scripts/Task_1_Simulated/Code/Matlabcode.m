    clc; clear; close all; 
    
    %% Loading ECG Data
    
    % Set working directory to ECG data folder
    cd('Data');
    
    % Load ECG record using WFDB toolbox 
    [signal, fs, siginfo] = rdsamp('100', [], 3600); % Load all channels, first 3600 samples (10s at 360Hz)
    
    % Signal is a matrix with 2 columns (2 ECG leads: MLII and V5)
    % We take column 1 (MLII lead) and transpose to a row vector
    ecgRaw = signal(:, 1)';
    
    % Total number of samples loaded
    bufferSize = length(ecgRaw);
    
    % Build time axis (converts sample numbers to seconds)
    t_axis = (0 : bufferSize-1) / fs;
    
    % Print summary to command window (comment out later!!!!!!!!!!!!!!!)
    fprintf('Sampling Rate : %d Hz\n', fs);
    fprintf('Samples loaded: %d (%.1f sec)\n', bufferSize, bufferSize/fs);               
    
    %% Signal Filtering
    
    % Filter Parameters
    % fs = 360 Hz (as defined above)
    notchFreq = 50; % Power line noise frequency in Pakistan
    
    % Filter 1: Notch filter (removes 50 Hz power line noise) 
    % iirnotch creates a very narrow notch exactly at 50 Hz
    % wo = normalized frequency (50 divided by Nyquist frequency fs/2)
    % bw = width of the notch (narrow so we only remove 50 Hz)
    wo = notchFreq / (fs/2);
    bw = wo / 35;
    [n_b, n_a] = iirnotch(wo, bw);
    % filtfilt applies filter TWICE (forward + backward) and zero phase distortion
    % meaning R-peak positions stay at their correct time locations
    ecgNotched  = filtfilt(n_b, n_a, ecgRaw);
    
    % Filter 2: Bandpass filter (0.5 Hz to 40 Hz) ---
    % High-pass at 0.5 Hz  → removes slow baseline wander
    % Low-pass  at 40 Hz   → removes high frequency muscle noise
    % butter(2, ...) = 2nd order Butterworth (smooth, no ripple)
    % [0.5 40]/(fs/2) = normalize cutoff frequencies to Nyquist
    [bp_b, bp_a]  = butter(2, [0.5 40] / (fs/2), 'bandpass');
    
    % Apply bandpass filter (again filtfilt for zero phase)
    ecgFiltered   = filtfilt(bp_b, bp_a, ecgNotched);
    
    %% R-Peak Detection
    
    % Normalize the filtered signal ---
    % Normalize to range [-1, 1] so our threshold works regardless
    % of the original amplitude scale of the signal
    ecgNorm = (ecgFiltered - mean(ecgFiltered)) / max(abs(ecgFiltered));
    
    % Set detection parameters
    % Threshold: R-peaks must be above 50% of the maximum value
    % This filters out smaller waves (P and T waves) which are lower
    threshold = 0.5 * max(ecgNorm);
    
    % Minimum distance between two R-peaks (for healthy heart 0.3 sec between beats)
    % 0.3 sec * 360 Hz = 108 samples minimum distance
    minPeakDist = round(0.3 * fs);
    
    % Detect R-peaks using findpeaks
    % MinPeakHeight : ignore anything below threshold
    % MinPeakDistance : ignore peaks too close together
    [peakVals, peakLocs] = findpeaks(ecgNorm, ...
        'MinPeakHeight',   threshold, ...
        'MinPeakDistance', minPeakDist);
    
    %% Heart Rate (BPM) and HRV Calculations
    
    % Step 1: Calculate RR Intervals 
    % RR interval = time between consecutive R-peaks
    RR_samples  = diff(peakLocs); % distance between peaks in samples
    RR_sec      = RR_samples / fs; % convert to seconds
    RR_ms       = RR_sec * 1000; % convert to milliseconds
    
    % Step 2: Heart Rate (BPM)
    % Formula: BPM = 60 / average RR interval in seconds
    BPM         = 60 / mean(RR_sec);
    
    % Step 3: HRV Metrics 
    
    % SDNN: Standard Deviation of RR intervals
    % Measures overall variability — higher = healthier autonomic system
    HRV_SDNN    = std(RR_ms);
    
    % RMSSD: Root Mean Square of Successive Differences
    % Measures short-term beat-to-beat variability
    HRV_RMSSD   = sqrt(mean(diff(RR_ms).^2));
    
    % pNN50: Percentage of successive RR differences greater than 50ms
    % High pNN50 = good parasympathetic (rest) nervous system activity
    pNN50_count = sum(abs(diff(RR_ms)) > 50);   % count pairs differing > 50ms
    HRV_pNN50   = (pNN50_count / length(diff(RR_ms))) * 100;  % convert to %
    
    %% ECG Dashboard 
    
    % Extract one complete PQRST cycle
    preWin  = round(0.3 * fs);
    postWin = round(0.5 * fs);
    
    if length(peakLocs) >= 4
        idx = peakLocs(3);
        if (idx - preWin >= 1) && (idx + postWin <= length(ecgNorm))
            oneCycle = ecgNorm(idx - preWin : idx + postWin);
            cycle_t  = (-preWin : postWin) / fs * 1000;   % ms
        end
    else
        oneCycle = [];
        warning('Not enough R-peaks to extract a single cycle.');
    end
    
    % Build Dashboard Figure 
    fig = figure('Name',        'ECG Dashboard', ...
                 'NumberTitle', 'off', ...
                 'Position',    [30 30 1100 600], ...
                 'Color',       [0.82 0.82 0.82]);   % grey outer background
    
    % PANEL 1: ECG signal in real time
    ax1 = subplot(2, 2, [1 2]);
    plot(t_axis, ecgRaw, ...
        'Color',    [0.75 0.75 0.75], ...
        'LineWidth', 0.8);
    hold on;
    
    plot(t_axis, ecgNorm, ...
        'Color',    [0.18 0.80 0.44], ...
        'LineWidth', 1.4);
    
    plot(peakLocs/fs, peakVals, ...
        'rv', ...
        'MarkerSize',      5, ...
        'MarkerFaceColor', 'r');
    
    legend('Raw ECG', 'Filtered ECG', 'R-peaks', ...
           'Location', 'northeast', ...
           'FontSize',  8);
    title('Full ECG Waveform', 'FontSize', 13);
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    grid on;
    hold off;
    
    
    % PANEL 2: ECG signal one cycle
    ax2 = subplot(2, 2, 3);
    if ~isempty(oneCycle)
        area(cycle_t, oneCycle, ...
            'FaceColor', [0.10 0.45 0.85], ...
            'FaceAlpha',  0.25, ...
            'EdgeColor', [0.30 0.65 1.0], ...
            'LineWidth',  1.8);
        hold on;
    
        % Vertical line at R-peak
        yl = ylim;
        line([0 0], yl, ...
            'LineStyle', '--', ...
            'Color',     'r', ...
            'LineWidth',  1.2);
        text(5, yl(2) * 0.90, 'R', ...
            'Color',      'r', ...
            'FontWeight', 'bold', ...
            'FontSize',    6);
    
        % Label PQRST waves
        [~, rIdx] = max(oneCycle);
        text(cycle_t(rIdx), oneCycle(rIdx) + 0.05, 'R', ...
            'Color',               'r', ...
            'FontWeight',          'bold', ...
            'HorizontalAlignment', 'center');
        text(-250, 0.3, 'P', 'Color', 'c', 'FontWeight', 'bold');
        text( 200, 0.2, 'T', 'Color', 'y', 'FontWeight', 'bold');
    
        xlabel('Time (ms)');
        ylabel('Amplitude');
        hold off;
    else
        text(0.3, 0.5, 'Insufficient peaks', 'Units', 'normalized');
    end
    title('Single ECG Cycle', 'FontSize', 12);
    grid on;
    
    % PANEL 3: Metrics display
    ax3 = subplot(2, 2, 4);
    
    set(ax3, 'Color', [1 1 1], ...
             'LineWidth', 1.5, ...
             'Box', 'on');
    axis off;
    
    % Heart Rate
    text(0.10, 0.72, 'Heart Rate =', ...
        'FontSize',   14, ...
        'FontWeight', 'bold', ...
        'Color',      [0 0 0], ...
        'Units',      'normalized');
    text(0.72, 0.72, sprintf('%.1f BPM', BPM), ...
        'FontSize',   14, ...
        'FontWeight', 'bold', ...
        'Color',      [0 0 0], ...
        'Units',      'normalized');
    
    % Divider line
    annotation('line', ...
        ax3.Position(1) + [0.02 0.96] * ax3.Position(3), ...
        [1 1] * (ax3.Position(2) + 0.52 * ax3.Position(4)), ...
        'Color',     [0.6 0.6 0.6], ...
        'LineWidth',  0.8);
    
    % HRV (SDNN)
    text(0.10, 0.40, 'HRV =', ...
        'FontSize',   14, ...
        'FontWeight', 'bold', ...
        'Color',      [0 0 0], ...
        'Units',      'normalized');
    text(0.48, 0.40, sprintf('%.1f ms', HRV_SDNN), ...
        'FontSize',   14, ...
        'FontWeight', 'bold', ...
        'Color',      [0 0 0], ...
        'Units',      'normalized');
    
    % Sub-metrics in smaller text
    text(0.10, 0.20, sprintf('RMSSD: %.2f ms     pNN50: %.1f%%', HRV_RMSSD, HRV_pNN50), ...
        'FontSize',  9, ...
        'Color',     [0.35 0.35 0.35], ...
        'Units',     'normalized');
    
    text(0.10, 0.08, sprintf('R-peaks: %d     Duration: %.0f sec', length(peakLocs), bufferSize/fs), ...
        'FontSize',  9, ...
        'Color',     [0.35 0.35 0.35], ...
        'Units',     'normalized');
    
    % Overall title
    sgtitle(sprintf('ECG Dashboard |  %s', ...
            datestr(now, 'dd-mmm-yyyy HH:MM')), ...
            'FontSize', 13, 'FontWeight', 'bold');

