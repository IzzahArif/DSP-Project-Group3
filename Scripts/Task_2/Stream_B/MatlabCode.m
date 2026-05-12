%B.1 — Artifact Dataset Selection:
cd('C:\Users\ramee\Desktop\University\Semester 6\DSP\Project\ECG_Data'); % Load clean ECG
[clean_signal, fs] = rdsamp('118', [], 3600);
ecg_clean_ref = clean_signal(:,1)';
[noise_signal, fs_n] = rdsamp('em', [], 3600);             % Load electrode motion noise
em_noise = noise_signal(:,1)';
minLen = min(length(ecg_clean_ref), length(em_noise));     % Match lengths
ecg_clean_ref = ecg_clean_ref(1:minLen);
em_noise      = em_noise(1:minLen);
% Add known motion artifacts to a clean signal at various SNR levels,
target_SNRs = [20, 10, 5];   % dB
ecg_contaminated = cell(1,3);
signal_power = mean(ecg_clean_ref.^2);
noise_power  = mean(em_noise.^2);
for k = 1:3
   scale = sqrt(signal_power / (noise_power * 10^(target_SNRs(k)/10)));
   ecg_contaminated{k} = ecg_clean_ref + scale * em_noise;
end
% Compute scale for SNR = 10 dB (k=2) explicitly
scale_10dB = sqrt(signal_power / (noise_power * 10^(10/10)));
ecg_artifact = ecg_contaminated{2};  % Work with SNR = 10 dB for the rest of Stream B
t_axis = (0:minLen-1)/fs;
% Plot
figure;
subplot(3,1,1); 
plot(t_axis, ecg_clean_ref, 'b');
title('Clean ECG Reference'); 
xlabel('Time (s)'); ylabel('mV'); grid on;
subplot(3,1,2); 
plot(t_axis, em_noise, 'r');
title('Electrode Motion Noise');
 xlabel('Time (s)'); ylabel('mV'); grid on;
subplot(3,1,3);
 plot(t_axis, ecg_artifact, 'k');
title('Contaminated ECG (SNR = 10 dB)');
 xlabel('Time (s)'); ylabel('mV'); grid on;
%  Document the artifact characteristics.
fprintf('=== Artifact Characteristics ===\n');
fprintf('Clean record      : MIT-BIH 118\n');
fprintf('Noise source      : MIT-BIH Noise Stress Test DB (em)\n');
fprintf('Sampling Rate     : %d Hz\n', fs);
fprintf('Signal Duration   : %.2f sec\n', minLen/fs);
fprintf('Artifact SNR used : 10 dB\n');
%B.2 — Demonstrate Fixed Filter Failure
% Apply Stream A filtering pipeline (HP -> Notch -> LP)
hp_cutoff  = 0.5;                                            % high pass 
[b_hp, a_hp] = butter(2, hp_cutoff/(fs/2), 'high');
notchFreq = 50;                                              % notch 
wo = notchFreq/(fs/2);  bw = wo/35;
[n_b, n_a] = iirnotch(wo, bw);
lp_cutoff = 40;                                              % low pass
[b_lp, a_lp] = butter(4, lp_cutoff/(fs/2), 'low');
% Apply fixed pipeline
ecg_fixed = filtfilt(b_hp, a_hp, ecg_artifact);
ecg_fixed = filtfilt(n_b,  n_a,  ecg_fixed);
ecg_fixed = filtfilt(b_lp, a_lp, ecg_fixed);
% Show failure
t_start = 5; t_end = 10;                          % adjust to where artifact is visible
idx_zoom = (t_axis >= t_start & t_axis <= t_end);
figure;
subplot(3,1,1); 
plot(t_axis(idx_zoom), ecg_clean_ref(idx_zoom), 'b');
title('Clean ECG (Reference)');
xlabel('Time (s)'); ylabel('mV'); grid on;
subplot(3,1,2); 
plot(t_axis(idx_zoom), ecg_artifact(idx_zoom), 'r');
title('Contaminated ECG');
xlabel('Time (s)'); ylabel('mV'); grid on;
subplot(3,1,3); 
plot(t_axis(idx_zoom), ecg_fixed(idx_zoom), 'm');
title('After Fixed Filter Pipeline — Ringing / Distortion Visible');
xlabel('Time (s)'); ylabel('mV'); grid on;
% SNR of fixed filter output
P_sig_fixed   = mean((ecg_fixed   - mean(ecg_fixed)).^2);
P_noise_fixed = mean((ecg_fixed   - ecg_clean_ref).^2);
SNR_fixed     = 10*log10(P_sig_fixed / P_noise_fixed);
fprintf('\n=== Fixed Filter SNR = %.2f dB ===\n', SNR_fixed);
%B.3 — Time-Frequency Analysis (STFT / Spectrogram):
% Compute the Short-Time Fourier Transform (STFT) of the artifact-contaminated signal.
win_sizes = [64, 256, 512];   % different window size for comparison
win_labels = {'64 samples (good time res)','256 samples (balanced)','512 samples (good freq res)'};
figure;
for w = 1:3
   win = hamming(win_sizes(w));
   noverlap = floor(win_sizes(w)/2);
   nfft = max(256, win_sizes(w));
   subplot(3,1,w);
   spectrogram(ecg_artifact, win, noverlap, nfft, fs, 'yaxis');
   title(['Spectrogram — Window = ' win_labels{w}]);
   ylim([0 80]);
   colorbar;
end
sgtitle('STFT Spectrogram: Motion Artifact Signal');
% Plot the working spectrogram (time vs. frequency vs. magnitude).
% time (x axis) , freq (y axis) , magnitude (colour)
figure;
spectrogram(ecg_artifact, hamming(256), 128, 512, fs, 'yaxis');
title('Spectrogram — Contaminated ECG (256-pt window)');
ylim([0 80]); colorbar;
% Identify how artifact energy appears differently from ECG energy in the time-frequency plane.
% ECG energy mainly appears below 40 Hz with regular heartbeat patterns, while motion 
% artifacts appear as irregular low-frequency bursts and EMG noise appears at higher 
% frequencies. 

% Experiment with different window sizes and discuss the time-frequency resolution 
% trade-off.
% Smaller windows give better time resolution but poorer frequency resolution, whereas 
% larger windows improve frequency resolution but blur sudden changes in time. 
%B.4 — Frame-Based / Segmented Processing:
% Divide the signal into overlapping frames (256 samples, 50% overlap)
% Apply windowing (Hamming)
% Extract features (energy, variance, kurtosis)
% Detect artifact frames using thresholding
frameLen  = 256;
hopSize   = 128;    % 50% overlap
win_func  = hamming(frameLen);
numFrames = floor((minLen - frameLen) / hopSize) + 1;
frame_energy   = zeros(1, numFrames);   % detects sudden motion spikes 
frame_variance = zeros(1, numFrames);   % very good for motion artifact detection (unstable signal) 
frame_kurtosis = zeros(1, numFrames);   % excellent for spike detection (motion artifacts are impulsive) 
frame_times    = zeros(1, numFrames);   % maps features back to time axis 
for i = 1:numFrames
   idx_start = (i-1)*hopSize + 1;
   idx_end   = idx_start + frameLen - 1;
   frame     = ecg_artifact(idx_start:idx_end) .* win_func';
   frame_energy(i)   = sum(frame.^2);
   frame_variance(i) = var(frame);
   frame_kurtosis(i) = kurtosis(frame);
   frame_times(i)    = (idx_start + frameLen/2 - 1) / fs;
end
% Plot frame-level statistics
figure;
subplot(3,1,1); plot(frame_times, frame_energy, 'b');
title('Frame Energy'); xlabel('Time (s)'); ylabel('Energy'); grid on;
subplot(3,1,2); plot(frame_times, frame_variance, 'r');
title('Frame Variance'); xlabel('Time (s)'); ylabel('Variance'); grid on;
subplot(3,1,3); plot(frame_times, frame_kurtosis, 'g');
title('Frame Kurtosis'); xlabel('Time (s)'); ylabel('Kurtosis'); grid on;
sgtitle('Frame-Based Signal Statistics');
%B.5 — Artifact Detection Algorithm:
% Implement an artifact detector that flags corrupted segments.
% Methods to consider:
% (a) Moving-average energy detector,
% (b) Variance-based threshold,
% (c) Wavelet-based detection.
% Mark artifact regions on your time-domain plot.
signalLen = length(ecg_artifact);
% (A) Moving-Average Energy Detector
energy_mean = mean(frame_energy);
energy_std  = std(frame_energy);
threshold_E = energy_mean + 2*energy_std;  % threshold
artifact_energy_frames = frame_energy > threshold_E;
artifact_mask_E = false(1, signalLen);
for i = 1:numFrames
   if artifact_energy_frames(i)
       idx_start = (i-1)*hopSize + 1;
       idx_end   = min(idx_start + frameLen - 1, signalLen);
       artifact_mask_E(idx_start:idx_end) = true;
   end
end
% (B) Variance-Based Detector
var_mean = mean(frame_variance);
var_std  = std(frame_variance);
threshold_V = var_mean + 2*var_std;  % treshold
artifact_var_frames = frame_variance > threshold_V;
artifact_mask_V = false(1, signalLen);
for i = 1:numFrames
   if artifact_var_frames(i)
       idx_start = (i-1)*hopSize + 1;
       idx_end   = min(idx_start + frameLen - 1, signalLen);
       artifact_mask_V(idx_start:idx_end) = true;
   end
end
% (C) Wavelet-Based Detection
[c,l] = wavedec(ecg_artifact, 5, 'db4');
% Use detail coefficients (high frequency components)
d3 = detcoef(c,l,3);
wave_energy = abs(d3).^2;
wave_mean = mean(wave_energy);
wave_std  = std(wave_energy);
threshold_W = wave_mean + 2*wave_std;  % threshold
artifact_wave_frames = wave_energy > threshold_W;
artifact_mask_W = false(1, signalLen);
stepW = floor(signalLen/length(wave_energy));
for i = 1:length(artifact_wave_frames)
   if artifact_wave_frames(i)
       idx_start = (i-1)*stepW + 1;
       idx_end   = min(i*stepW, signalLen);
       artifact_mask_W(idx_start:idx_end) = true;
   end
end
% FINAL COMBINED DETECTION
artifact_mask_final = artifact_mask_E | artifact_mask_V | artifact_mask_W;
% PLOT
figure;
plot(t_axis, ecg_artifact, 'b'); hold on;
% Highlight detected regions
plot(t_axis(artifact_mask_E), ecg_artifact(artifact_mask_E), 'r.', 'MarkerSize', 2);
plot(t_axis(artifact_mask_V), ecg_artifact(artifact_mask_V), 'g.', 'MarkerSize', 2);
plot(t_axis(artifact_mask_W), ecg_artifact(artifact_mask_W), 'k.', 'MarkerSize', 2);
title('B.5 Artifact Detection (Energy / Variance / Wavelet)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
legend('ECG','Energy','Variance','Wavelet');
grid on;
% Print Results
fprintf('\n========== ARTIFACT DETECTION ==========\n');
fprintf('Energy Threshold     : %.4f\n', threshold_E);
fprintf('Variance Threshold   : %.4f\n', threshold_V);
fprintf('Wavelet Threshold    : %.4f\n', threshold_W);
fprintf('Energy Artifacts     : %d frames\n', sum(artifact_energy_frames));
fprintf('Variance Artifacts   : %d frames\n', sum(artifact_var_frames));
fprintf('Wavelet Artifacts    : %d frames\n', sum(artifact_wave_frames));
fprintf('Total Artifact Cover : %.2f %%\n', 100*mean(artifact_mask_final));
%B.6 — Artifact Suppression / Removal:
% For detected artifact segments, apply one or more of:
% (a) Interpolation — replace corrupted segment with interpolated values from clean neighbors,
% (b) Adaptive filtering — use LMS or RLS adaptive filter with a reference noise channel,
% (c) Median filtering — for spike removal,
% (d) Wavelet denoising — decompose, threshold artifact coefficients, reconstruct.
% Show before/after results for each method attempted.
% (a) Interpolation
ecg_interp = ecg_artifact;
artifact_regions = bwlabel(artifact_mask_final);   % label connected regions
num_regions = max(artifact_regions);
for r = 1:num_regions
   region_idx = find(artifact_regions == r);
   left_idx   = max(1,      region_idx(1) - 1);
   right_idx  = min(minLen, region_idx(end) + 1);
   left_val   = ecg_artifact(left_idx);
   right_val  = ecg_artifact(right_idx);
   interp_vals = linspace(left_val, right_val, length(region_idx));
   ecg_interp(region_idx) = interp_vals;
end
% (b) Adaptive LMS Filter
mu         = 0.001;       % LMS step size
filter_len = 32;         % adaptive filter length
ecg_lms    = zeros(1, minLen);
w_lms      = zeros(1, filter_len);   % filter weights
ref = scale_10dB * em_noise;         % reference noise (from B.1 scaling)
for n = filter_len:minLen
   x_vec       = ref(n:-1:n-filter_len+1);    % reference window
   y_hat       = w_lms * x_vec';              % filter output (noise estimate)
   e           = ecg_artifact(n) - y_hat;     % error = cleaned sample
   w_lms       = w_lms + 2*mu*e*x_vec;        % LMS weight update
   ecg_lms(n)  = e;
end
% (c) Median filtering
ecg_median = medfilt1(ecg_artifact, 5);   % 5-point median filter
% (d) Wavelet Denoising
wname  = 'db4';
level  = 5;
[C, L] = wavedec(ecg_artifact, level, wname);
sigma = median(abs(C)) / 0.6745;   % Noise estimation (robust)
thr = sigma * sqrt(2 * log(minLen));   % universal threshold
C_thresh = wthresh(C, 's', thr);       % soft thresholding
ecg_wavelet = waverec(C_thresh, L, wname);
% PLOTS
figure;
subplot(5,1,1);
plot(t_axis, ecg_artifact, 'b');
title('Original Artifact ECG');
ylabel('mV'); grid on;
subplot(5,1,2);
plot(t_axis, ecg_interp, 'g');
title('(a) Interpolation Method');
ylabel('mV'); grid on;
subplot(5,1,3);
plot(t_axis, ecg_lms, 'm');
title('(b) LMS Adaptive Filtering');
ylabel('mV'); grid on;
subplot(5,1,4);
plot(t_axis, ecg_median, 'c');
title('(c) Median Filtering');
ylabel('mV'); grid on;
subplot(5,1,5);
plot(t_axis, ecg_wavelet, 'r');
title('(d) Wavelet Denoising');
ylabel('mV');
xlabel('Time (s)');
grid on;
sgtitle('Artifact Removal — Method Comparison');
% Display Results
if exist('ecg_clean_ref','var')
   rmse_interp = sqrt(mean((ecg_clean_ref - ecg_interp).^2));
   rmse_lms    = sqrt(mean((ecg_clean_ref - ecg_lms).^2));
   rmse_med    = sqrt(mean((ecg_clean_ref - ecg_median).^2));
   rmse_wave   = sqrt(mean((ecg_clean_ref - ecg_wavelet).^2));
   fprintf('\n========== PERFORMANCE ==========\n');
   fprintf('RMSE Interpolation : %.4f\n', rmse_interp);
   fprintf('RMSE LMS           : %.4f\n', rmse_lms);
   fprintf('RMSE Median        : %.4f\n', rmse_med);
   fprintf('RMSE Wavelet       : %.4f\n', rmse_wave);
end
%B.7 — Performance Comparison:
% Compare at least two artifact removal methods quantitatively.
% Use metrics: SNR improvement, waveform correlation with clean reference, RMS error.
% Helper: compute SNR and RMSE vs clean reference
compute_snr  = @(sig, ref) 10*log10(mean(ref.^2) / mean((sig-ref).^2));
compute_rmse = @(sig, ref) sqrt(mean((sig-ref).^2));
compute_corr = @(sig, ref) corr(sig(:), ref(:));
methods = {'Contaminated','Fixed Filter','Interpolation','LMS Adaptive','Median','Wavelet'};
signals = {ecg_artifact, ecg_fixed, ecg_interp, ecg_lms, ecg_median, ecg_wavelet};
SNR_vals  = zeros(1,6);
RMSE_vals = zeros(1,6);
CORR_vals = zeros(1,6);
for m = 1:6
    SNR_vals(m)  = compute_snr( signals{m}, ecg_clean_ref);
    RMSE_vals(m) = compute_rmse(signals{m}, ecg_clean_ref);
    CORR_vals(m) = compute_corr(signals{m}, ecg_clean_ref);
end
fprintf('\n=== B.7 Performance Comparison Table ===\n');
fprintf('%-20s %10s %10s %10s\n', 'Method', 'SNR (dB)', 'RMSE', 'Correlation');
fprintf('%s\n', repmat('-',1,55));
for m = 1:6
    fprintf('%-20s %10.2f %10.4f %10.4f\n', ...
        methods{m}, SNR_vals(m), RMSE_vals(m), CORR_vals(m));
end
% Create a comparison table.
figure;
subplot(1,3,1); bar(SNR_vals);
set(gca,'xticklabels', methods,'XTickLabelRotation', 30);
title('SNR (dB)'); ylabel('dB'); grid on;
subplot(1,3,2); bar(RMSE_vals);
set(gca,'xticklabels', methods,'XTickLabelRotation', 30);
title('RMSE'); ylabel('mV'); grid on;
subplot(1,3,3); bar(CORR_vals);
set(gca,'xticklabels', methods,'XTickLabelRotation', 30);
title('Correlation with Clean'); ylabel('r'); ylim([0 1]); grid on;
sgtitle('B.7 — Quantitative Comparison of Artifact Removal Methods');


