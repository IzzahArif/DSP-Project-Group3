% A.1 — Signal Loading & Time-Domain Inspection:
clear; clc; fclose('all'); close all;
cd('C:\Users\izzah\Desktop\Project\ECG_Data');     % Record ID: 100
% Manual decoding of MIT-BIH Format 212
fid = fopen('100.dat', 'r');
if fid == -1
   error('File 100.dat not found! Ensure it is in the specified folder.');
end
A = fread(fid, [3, inf], 'uint8')';
fclose(fid);
% Decoding bit-packed samples
M = A(:, 2);
sig1 = bitshift(bitand(M, 15), 8) + A(:, 1);
sig1(sig1 > 2047) = sig1(sig1 > 2047) - 4096; % Signed 12-bit conversion
fs = 360;          
gain = 200;        
base = 1024;       
ecgRaw = (sig1' - base) / gain;               % Convert to mV
bufferSize = length(ecgRaw);
t_axis = (0:bufferSize-1)/fs;
% Plot Raw signal vs Time
figure('Name', 'A.1 Raw ECG Signal');
plot(t_axis, ecgRaw, 'b');
title('Raw ECG Signal (Record 100)');
xlabel('Time (s)'); ylabel('Amplitude (mV)'); grid on;
xlim([0 10]);

% Identify visible noise artifacts.
% Baseline wander (slow movement of ECG baseline), power-line interference (50 Hz sinusoidal noise),
% EMG noise (muscle movement), random spikes/artifacts (motion or electrodes)
% Document the sampling rate, duration, and record ID.
fprintf('Record ID        : 100\n');
fprintf('Sampling Rate    : %d Hz\n', fs);
fprintf('Signal Duration  : %.2f sec\n', bufferSize/fs);
% A.2 — FFT & Frequency-Domain Analysis:
% FFT of the raw ECG signal
N = length(ecgRaw);
X = fft(ecgRaw);
X_mag = abs(X)/N;                         % normalized magnitude
% single-sided magnitude spectrum 
f = (0:N-1)*(fs/N);
half_N = floor(N/2);
f_half = f(1:half_N);
X_half = 2*X_mag(1:half_N);
figure('Name', 'A.2 Magnitude Spectrum');
plot(f_half, X_half, 'LineWidth', 1.5);
title('Single-Sided ECG Spectrum');
xlabel('Frequency (Hz)'); ylabel('|X(f)|'); xlim([0 150]); grid on;
ylim([0 0.05]) 
hold on;
xline(0.5,'g--','0.5 Hz');              % ECG energy band
xline(40,'g--','40 Hz');
xline(60,'r--','60 Hz Noise');          % Power-line frequency
xline(100,'y--','100 Hz');                % High frequency noise
hold off;
% Identify peaks:
[pks, locs] = findpeaks(X_half, f_half,'MinPeakHeight',0.02);
hold on;
plot(locs,pks,'ro');
for i=1:length(pks)
   text(locs(i),pks(i),sprintf('%.1f Hz',locs(i)));
end
hold off;
% A.3 — Noise Power Estimation:
% power in each noise band vs. signal band.
signal_idx = (f_half >= 0.5 & f_half <= 40); 
P_signal = sum(X_half(signal_idx).^2);      % signal power (0.5–40 Hz)
bw_idx = (f_half < 0.5);     
P_bw = sum(X_half(bw_idx).^2);              % baseline wander power (< 0.5 Hz)
pl_idx = (f_half >= 59 & f_half <= 61);    
P_pl = sum(X_half(pl_idx).^2);              % power line noise power (59 - 61 Hz)
emg_idx = (f_half > 100);
P_emg = sum(X_half(emg_idx).^2);            % emg noise power (>100) 
P_noise = P_bw + P_pl + P_emg;               % TOTAL NOISE POWER

% Compute SNR before filtering.
SNR_before = 10*log10(P_signal/P_noise);
fprintf('\n========== SNR BEFORE FILTERING ==========\n');
fprintf('Signal Power = %.4f\n', P_signal);
fprintf('Noise Power  = %.4f\n', P_noise);
fprintf('SNR before   = %.2f dB\n', SNR_before);
% A.4 — Filter Design for Baseline Wander Removal:
fir_order_hp = 80; hp_order = 4; 
[b_hp, a_hp] = butter(hp_order, hp_cutoff/(fs/2), 'high'); % high-Pass IIR filter
b_fir = fir1(fir_order_hp, hp_cutoff/(fs/2), 'high');      % high-Pass FIR filter
figure;
freqz(b_hp, a_hp, 1024, fs);                                 % IIR filter freq response
ecg_hp = filtfilt(b_hp, a_hp, ecgRaw); 
figure;
freqz(b_fir,1,1024,fs);                                      % FIR filter freq response
ecg_fir = filtfilt(b_fir,1,ecgRaw);

% FIR filters have linear phase, so all frequencies are delayed equally and ECG shape is preserved.
% IIR filters have non-linear phase, so different frequencies are delayed differently, which may slightly 
% distort the waveform but they are more efficient. 

% A.5 — Filter Design for Power-Line Interference:
notchFreq = 60;                                    % notch filter at 60 Hz 
wo = notchFreq/(fs/2);
bw = wo/35;                                        % 35 is tuning factor
[n_b, n_a] = iirnotch(wo, bw);                     % IIR notch filter 
% magnitude response showing the notch.
figure;
freqz(n_b,n_a,1024,fs);
% Apply the filter
ecg_notch = filtfilt(n_b,n_a,ecg_hp); 
% Verify removal in the frequency spectrum.
figure;
subplot(2,1,1);
plot(f_half,X_half);
title('Before Notch');
subplot(2,1,2);
X2 = abs(fft(ecg_notch))/N;
plot(f_half,2*X2(1:half_N));
title('After Notch');

% A.6 — Filter Design for EMG Noise:
% Design a low-pass filter with cutoff ~40–100 Hz 
% ECG lies below 40 Hz while EMG noise is mainly above it. 
lp_cutoff = 40;              
[b_lp,a_lp] = butter(4, lp_cutoff/(fs/2), 'low');           % IIR Filter with order 4 
figure;
freqz(b_lp,a_lp,1024,fs);                                   % freq response
ecg_lp = filtfilt(b_lp,a_lp,ecg_notch);                     % apply filter
fir_order_lp = 80;                                          % FIR filter with order 80 
b_fir_lp = fir1(fir_order_lp, lp_cutoff/(fs/2), 'low');
figure;
freqz(b_fir_lp,1,1024,fs);
% Comparison FIR vs. IIR. 
figure;
subplot(2,1,1);
plot(f_half,X_half);
title('Before LP Filter');
xlim([0 150]);
subplot(2,1,2);
X_clean = abs(fft(ecg_lp))/N;
plot(f_half,2*X_clean(1:half_N));
title('After LP Filter');
xlim([0 150]);

% A.7 — Combined Filtering Pipeline:
ecg_clean = filtfilt(b_hp,a_hp,ecgRaw);
ecg_clean = filtfilt(n_b,n_a,ecg_clean);
ecg_clean = filtfilt(b_lp,a_lp,ecg_clean);
ecg_clean_fir = filtfilt(b_fir,1,ecgRaw);
ecg_clean_fir = filtfilt(n_b,n_a,ecg_clean_fir); 
ecg_clean_fir = filtfilt(b_fir_lp,1,ecg_clean_fir);
% Plot:
% (a) raw vs. cleaned signal in time domain,
figure;
subplot(2,1,1); plot(t_axis, ecgRaw); title('Raw ECG');
subplot(2,1,2); plot(t_axis, ecg_clean); title('Cleaned ECG');
% (b) raw vs. cleaned spectrum in frequency domain.
X_raw = abs(fft(ecgRaw))/N;
X_clean = abs(fft(ecg_clean))/N;
figure;
subplot(2,1,1); 
plot(f_half,2*X_raw(1:half_N)); title('Raw Spectrum'); xlim([0 150]);
subplot(2,1,2); 
plot(f_half,2*X_clean(1:half_N)); title('Filtered Spectrum'); xlim([0 150]);
% Comparison between FIR and IIR
figure;
plot(t_axis, ecg_clean, 'b'); hold on;
plot(t_axis, ecg_clean_fir, 'r');
legend('IIR Cleaned','FIR Cleaned');
title('FIR vs IIR ECG Comparison');
grid on;
% Compute SNR after filtering. Show SNR improvement.
Xf = abs(fft(ecg_clean))/N;
Xf_half = 2*Xf(1:half_N);
P_signal_after = sum(Xf_half(signal_idx).^2);
P_noise_after = sum(Xf_half(bw_idx).^2) + ...
   sum(Xf_half(pl_idx).^2) + ...
   sum(Xf_half(emg_idx).^2);
SNR_after = 10*log10(P_signal_after/P_noise_after);
fprintf('\n========== SNR AFTER FILTERING ==========\n');
fprintf('SNR After Filtering = %.2f dB\n', SNR_after);
fprintf('SNR Improvement = %.2f dB\n', SNR_after - SNR_before);

% A.8 — Group Delay Analysis:
% Compute and plot the group delay of each filter.
figure;                                                % high pass filter
grpdelay(b_hp,a_hp,1024,fs);
title('Group Delay — High Pass Filter');
figure;                                                % notch filter
grpdelay(n_b,n_a,1024,fs);
title('Group Delay — Notch Filter');
figure;                                                % low pass filter
grpdelay(b_lp,a_lp,1024,fs);
title('Group Delay — Low Pass Filter');
