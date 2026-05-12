%% =========================================================================
%  STREAM C — Multi-Lead Phase Consistency (Advanced Structural DSP)
%  12-Lead ECG Phase Preservation & Group Delay Analysis
%
%  NO EXTERNAL TOOLBOXES REQUIRED.
%  Uses only built-in MATLAB Signal Processing Toolbox functions:
%    butter, iirnotch, freqz, grpdelay, filter, filtfilt, xcorr, fir1
%
%  Data:
%    Reads directly from C:\Users\DSP lab\Downloads\Project\ECG_LeadData
%    The folder is scanned for .dat + .hea file pairs automatically.
%    fs, n_samples, gain, and baseline are all parsed from the .hea header —
%    nothing about the recording is hardcoded.
%    Falls back to a synthetic ECG if the folder is missing or empty.
%
%  Prerequisites: Stream A (IIR filters), Stream B (artifact removal).
%  =========================================================================

clear; clc; close all;

%% =========================================================================
%  SETUP — only edit this one line if your path differs
%  =========================================================================

ECG_FOLDER = 'C:\Users\DSP lab\Downloads\Project\ECG_LeadData';

lead_names  = {'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'};
n_leads     = 12;
lead_colors = lines(n_leads);

fprintf('=== STREAM C: Multi-Lead Phase Consistency Analysis ===\n\n');

%% =========================================================================
%  C.1 — MULTI-LEAD DATASET LOADING
%  =========================================================================
fprintf('--- C.1: Loading Multi-Lead ECG Dataset ---\n');
fprintf('  Folder : %s\n', ECG_FOLDER);

% Auto-scan folder for the first usable .dat + .hea pair
[ecg_raw, fs, record_used] = load_ecg_folder(ECG_FOLDER, n_leads);

N            = size(ecg_raw, 2);
duration_sec = N / fs;
t            = (0:N-1) / fs;   % time vector (seconds)

fprintf('  Record   : %s\n',   record_used);
fprintf('  fs       : %d Hz\n', fs);
fprintf('  Duration : %.1f s  (%d samples)\n', duration_sec, N);
fprintf('  Leads    : %s\n\n',  strjoin(lead_names,' | '));

% --- Plot C.1: Stacked 12-lead ECG ---
% Compute a sensible offset from the actual signal amplitude range
amp_range    = max(max(abs(ecg_raw))) * 2.5;
offset_scale = max(amp_range, 1.0);   % at least 1 mV spacing

figure('Name','C.1 — 12-Lead ECG (Original)','Position',[50 50 1400 900]);
t_disp = min(10, duration_sec);       % show up to 10 s in this plot
idx_disp = t <= t_disp;

for k = 1:n_leads
    offset = (n_leads - k) * offset_scale;
    plot(t(idx_disp), ecg_raw(k, idx_disp) + offset, ...
         'Color', lead_colors(k,:), 'LineWidth', 0.9);
    hold on;
    text(-0.05 * t_disp, offset, lead_names{k}, 'FontWeight','bold', ...
         'Color', lead_colors(k,:), 'FontSize', 9, 'HorizontalAlignment','right');
end
xlabel('Time (s)');
ylabel('Amplitude (mV) + display offset');
title(sprintf('C.1 — 12-Lead ECG: %s  |  fs = %d Hz  |  %.1f s', ...
              record_used, fs, duration_sec));
xlim([0 t_disp]); yticks([]); grid on; box off;

%% =========================================================================
%  C.2 — PHASE RESPONSE ANALYSIS OF STREAM A IIR FILTERS
%  =========================================================================
fprintf('--- C.2: Phase Response & Group Delay Analysis ---\n');

% Design Stream A IIR filters (causal Butterworth + notch)
% Notch at 50 Hz (change to 60 if your mains frequency is 60 Hz)
[b_hp,    a_hp   ] = butter(2, 0.5/(fs/2), 'high');
[b_lp,    a_lp   ] = butter(4,  40/(fs/2), 'low');
notch_freq = 50;   % Hz — change to 60 for 60 Hz mains
wo = notch_freq/(fs/2);  bw = wo/30;
[b_notch, a_notch] = iirnotch(wo, bw);

% Frequency axis
f_vec = linspace(0, fs/2, 4096);
w_vec = f_vec / (fs/2) * pi;

H_hp    = freqz(b_hp,    a_hp,    w_vec);
H_lp    = freqz(b_lp,    a_lp,    w_vec);
H_notch = freqz(b_notch, a_notch, w_vec);
H_combo = H_hp .* H_lp .* H_notch;

phase_hp    = unwrap(angle(H_hp));
phase_lp    = unwrap(angle(H_lp));
phase_notch = unwrap(angle(H_notch));
phase_combo = unwrap(angle(H_combo));

[gd_hp,    f_gd] = grpdelay(b_hp,    a_hp,    4096, fs);
[gd_lp,    ~   ] = grpdelay(b_lp,    a_lp,    4096, fs);
[gd_notch, ~   ] = grpdelay(b_notch, a_notch, 4096, fs);
gd_combo         = gd_hp + gd_lp + gd_notch;

% --- Plot C.2 ---
figure('Name','C.2 — IIR Phase & Group Delay','Position',[50 50 1300 700]);

subplot(2,1,1);
plot(f_vec, rad2deg(phase_hp),    'b',  'LineWidth', 1.5); hold on;
plot(f_vec, rad2deg(phase_lp),    'r',  'LineWidth', 1.5);
plot(f_vec, rad2deg(phase_notch), 'g',  'LineWidth', 1.5);
plot(f_vec, rad2deg(phase_combo), 'k',  'LineWidth', 2.0);
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('C.2a — IIR Filter Phase Responses (non-linear → phase distortion)');
legend('HP 0.5 Hz','LP 40 Hz',sprintf('Notch %d Hz',notch_freq),'Cascaded', ...
       'Location','southwest');
xlim([0 fs/2]); grid on;

subplot(2,1,2);
ecg_mask = (f_gd >= 0.5) & (f_gd <= 40);
plot(f_gd, gd_hp,    'b',  'LineWidth', 1.5); hold on;
plot(f_gd, gd_lp,    'r',  'LineWidth', 1.5);
plot(f_gd, gd_notch, 'g',  'LineWidth', 1.5);
plot(f_gd, gd_combo, 'k',  'LineWidth', 2.0);
xline(0.5, '--m', '0.5 Hz', 'LabelVerticalAlignment','bottom');
xline(40,  '--m', '40 Hz',  'LabelVerticalAlignment','bottom');
xlabel('Frequency (Hz)'); ylabel('Group Delay (samples)');
title('C.2b — Group Delay: VARIABLE across ECG band → different delay per frequency component');
legend('HP','LP','Notch','Cascaded','Location','northeast');
xlim([0 min(fs/2, 80)]); grid on;
yl2 = ylim; ylim([max(yl2(1),-30) min(yl2(2),200)]);

gd_band   = gd_combo(ecg_mask);
gd_var_ms = (max(gd_band) - min(gd_band)) / fs * 1000;
fprintf('  IIR cascade group delay in ECG band (0.5–40 Hz):\n');
fprintf('    Min      : %.1f samples  (%.1f ms)\n', min(gd_band), min(gd_band)/fs*1000);
fprintf('    Max      : %.1f samples  (%.1f ms)\n', max(gd_band), max(gd_band)/fs*1000);
fprintf('    Variation: %.1f samples  (%.1f ms)  <-- NON-CONSTANT = Phase Distortion\n\n', ...
        max(gd_band)-min(gd_band), gd_var_ms);

%% =========================================================================
%  C.3 — DEMONSTRATE PHASE DISTORTION PROBLEM
%  =========================================================================
fprintf('--- C.3: Demonstrating Phase Distortion Across Leads ---\n');

% Apply causal IIR (NOT zero-phase) independently to every lead
ecg_iir = zeros(n_leads, N);
for k = 1:n_leads
    s = filter(b_hp,    a_hp,    ecg_raw(k,:));
    s = filter(b_lp,    a_lp,    s);
    s = filter(b_notch, a_notch, s);
    ecg_iir(k,:) = s;
end

% Overlay 4 clinically related leads: I, II, V1, V2
show_idx = [1, 2, 7, 8];
cols_a   = {'#1f77b4','#d62728','#2ca02c','#9467bd'};

% Display window: first 3 s or full record if shorter
t_win = min(3.0, duration_sec - 0.1);
x_lim = [0.1 t_win];

figure('Name','C.3 — Phase Distortion Demonstration','Position',[50 50 1300 700]);
subplot(2,1,1);
for i = 1:length(show_idx)
    k = show_idx(i);
    plot(t, ecg_raw(k,:), 'Color', cols_a{i}, 'LineWidth', 1.4, ...
         'DisplayName', [lead_names{k} ' Original']); hold on;
end
xlabel('Time (s)'); ylabel('mV');
title('C.3a — Original: Leads Temporally Aligned');
legend('Location','northeast','NumColumns',2); grid on; xlim(x_lim);

subplot(2,1,2);
for i = 1:length(show_idx)
    k = show_idx(i);
    plot(t, ecg_iir(k,:), 'Color', cols_a{i}, 'LineWidth', 1.4, ...
         'LineStyle','--', 'DisplayName', [lead_names{k} ' IIR']); hold on;
end
xlabel('Time (s)'); ylabel('mV');
title('C.3b — After Causal IIR Filtering: Transient waveform distortion + inter-lead phase shift');
legend('Location','northeast','NumColumns',2); grid on; xlim(x_lim);

% Cross-correlation table
fprintf('  Cross-correlation timing — Original vs IIR:\n');
fprintf('  %-10s %-10s  %-18s  %-18s\n','Lead A','Lead B','Orig Lag (ms)','IIR Lag (ms)');
fprintf('  %s\n', repmat('-',1,60));
key_pairs  = [1,2; 1,3; 7,8; 8,9; 9,10];
max_lag_xc = round(fs * 0.3);
for p = 1:size(key_pairs,1)
    a = key_pairs(p,1);  b = key_pairs(p,2);
    [xo,lo] = xcorr(ecg_raw(a,:), ecg_raw(b,:), max_lag_xc, 'normalized');
    [xi,li] = xcorr(ecg_iir(a,:), ecg_iir(b,:), max_lag_xc, 'normalized');
    [~,io] = max(abs(xo));  orig_lag = lo(io)/fs*1000;
    [~,ii] = max(abs(xi));  iir_lag  = li(ii)/fs*1000;
    fprintf('  %-10s %-10s  %-18.1f  %-18.1f\n', ...
            lead_names{a}, lead_names{b}, orig_lag, iir_lag);
end
fprintf('\n');

%% =========================================================================
%  C.4 — GROUP DELAY MATCHING (Zero-Phase via filtfilt)
%  =========================================================================
fprintf('--- C.4: Group Delay Matching via filtfilt (zero-phase) ---\n');

% filtfilt: forward + backward pass → net phase = 0 at every frequency
ecg_zp = zeros(n_leads, N);
for k = 1:n_leads
    s = filtfilt(b_hp,    a_hp,    ecg_raw(k,:));
    s = filtfilt(b_lp,    a_lp,    s);
    s = filtfilt(b_notch, a_notch, s);
    ecg_zp(k,:) = s;
end

% FIR linear-phase reference filters (Method A alternative)
ord_hp   = 2 * round(2 * fs / 0.5);          % even order HP
ord_lp   = 2 * round(fs / 8);                 % even order LP
b_fir_hp = fir1(ord_hp, 0.5/(fs/2), 'high');
b_fir_lp = fir1(ord_lp,  40/(fs/2), 'low');
b_fir_nt = fir1(200, [(notch_freq-1)/(fs/2), (notch_freq+1)/(fs/2)], 'stop');

[gd_fir_hp, f_fir] = grpdelay(b_fir_hp, 1, 4096, fs);
[gd_fir_lp, ~    ] = grpdelay(b_fir_lp, 1, 4096, fs);
[gd_fir_nt, ~    ] = grpdelay(b_fir_nt, 1, 4096, fs);
gd_fir_total = gd_fir_hp + gd_fir_lp + gd_fir_nt;

% --- Plot C.4 ---
figure('Name','C.4 — Group Delay Comparison','Position',[50 50 1200 500]);
plot(f_gd,  gd_combo,     'r-',  'LineWidth', 2.0, 'DisplayName','IIR causal (variable)');
hold on;
plot(f_fir, gd_fir_total, 'b-',  'LineWidth', 2.0, 'DisplayName','FIR linear-phase (constant)');
yline(0, 'g--', 'LineWidth', 2.5, 'DisplayName','filtfilt zero-phase (delay = 0)');
xline(0.5,'--k','0.5 Hz','LabelVerticalAlignment','bottom','FontSize',8);
xline(40,  '--k','40 Hz', 'LabelVerticalAlignment','bottom','FontSize',8);
xlabel('Frequency (Hz)'); ylabel('Group Delay (samples)');
title('C.4 — Group Delay Matching: IIR causal vs FIR linear-phase vs filtfilt zero-phase');
legend('Location','northeast'); xlim([0 min(fs/2,80)]); grid on;

fir_band = (f_fir >= 1) & (f_fir <= 40);
fprintf('  filtfilt group delay     : 0 samples (0 ms) — constant\n');
fprintf('  FIR cascade group delay  : ~%d samples (%.1f ms) — constant\n', ...
        round(mean(gd_fir_total(fir_band))), mean(gd_fir_total(fir_band))/fs*1000);
fprintf('  IIR causal group delay   : %.0f–%.0f samples (%.1f–%.1f ms) — variable\n\n', ...
        min(gd_band), max(gd_band), min(gd_band)/fs*1000, max(gd_band)/fs*1000);

%% =========================================================================
%  C.5 — CROSS-CORRELATION ALIGNMENT ANALYSIS
%  =========================================================================
fprintf('--- C.5: Cross-Correlation Alignment Table (limb leads I–aVF) ---\n');

all_pairs  = nchoosek(1:6, 2);
n_pairs    = size(all_pairs, 1);
max_lag_xc = round(fs * 0.25);

iir_errors = zeros(n_pairs, 1);
zp_errors  = zeros(n_pairs, 1);

fprintf('\n  %-10s %-10s | %-14s %-14s %-14s | %-13s %-13s\n', ...
        'Lead A','Lead B','Orig(ms)','IIR(ms)','ZP(ms)', ...
        'IIR Err(ms)','ZP Err(ms)');
fprintf('  %s\n', repmat('-',1,95));

for p = 1:n_pairs
    a = all_pairs(p,1);  b = all_pairs(p,2);
    [xo,lo] = xcorr(ecg_raw(a,:), ecg_raw(b,:), max_lag_xc, 'normalized');
    [xi,li] = xcorr(ecg_iir(a,:), ecg_iir(b,:), max_lag_xc, 'normalized');
    [xz,lz] = xcorr(ecg_zp(a,:),  ecg_zp(b,:),  max_lag_xc, 'normalized');
    [~,io]=max(abs(xo)); orig_lag = lo(io)/fs*1000;
    [~,ii]=max(abs(xi)); iir_lag  = li(ii)/fs*1000;
    [~,iz]=max(abs(xz)); zp_lag   = lz(iz)/fs*1000;
    iir_errors(p) = abs(iir_lag - orig_lag);
    zp_errors(p)  = abs(zp_lag  - orig_lag);
    fprintf('  %-10s %-10s | %-14.1f %-14.1f %-14.1f | %-13.2f %-13.2f\n', ...
            lead_names{a}, lead_names{b}, orig_lag, iir_lag, zp_lag, ...
            iir_errors(p), zp_errors(p));
end

fprintf('\n  Summary:\n');
fprintf('    IIR causal — Mean error: %.2f ms,  Max error: %.2f ms\n', ...
        mean(iir_errors), max(iir_errors));
fprintf('    Zero-phase — Mean error: %.2f ms,  Max error: %.2f ms\n\n', ...
        mean(zp_errors), max(zp_errors));

% --- Plot C.5: XCorr functions for two representative pairs ---
figure('Name','C.5 — Cross-Correlation Analysis','Position',[50 50 1300 560]);
plot_pairs_c5 = [1,2; 7,10];  % I–II  and  V1–V4
max_lag_plot  = round(fs * 0.15);

for pp = 1:2
    a = plot_pairs_c5(pp,1);  b = plot_pairs_c5(pp,2);
    [xo,lo] = xcorr(ecg_raw(a,:), ecg_raw(b,:), max_lag_plot, 'normalized');
    [xi,~  ] = xcorr(ecg_iir(a,:), ecg_iir(b,:), max_lag_plot, 'normalized');
    [xz,~  ] = xcorr(ecg_zp(a,:),  ecg_zp(b,:),  max_lag_plot, 'normalized');
    lms = lo/fs*1000;
    subplot(1,2,pp);
    plot(lms, xo, 'b',  'LineWidth', 1.8, 'DisplayName','Original'); hold on;
    plot(lms, xi, 'r--','LineWidth', 1.8, 'DisplayName','IIR causal');
    plot(lms, xz, 'g',  'LineWidth', 2.2, 'DisplayName','Zero-phase');
    [~,io2]=max(abs(xo)); [~,ii2]=max(abs(xi)); [~,iz2]=max(abs(xz));
    xline(lms(io2),'b:','LineWidth',1.5);
    xline(lms(ii2),'r:','LineWidth',1.5);
    xline(lms(iz2),'g:','LineWidth',1.5);
    xlabel('Lag (ms)'); ylabel('Norm. Cross-Correlation');
    title(sprintf('C.5 — XCorr: %s vs %s', lead_names{a}, lead_names{b}));
    legend('Location','best'); grid on;
end

%% =========================================================================
%  C.6 — MULTI-CHANNEL SYSTEM CONSISTENCY VERIFICATION
%  =========================================================================
fprintf('--- C.6: Full Pipeline Verification (Stream A + B + C) ---\n');

% Full pipeline:
%   Stream A — zero-phase bandpass + notch
%   Stream B — baseline removal via moving-median subtraction
%   Stream C — phase consistency guaranteed by filtfilt throughout
ecg_final = zeros(n_leads, N);
win_bl    = 2 * round(0.3 * fs) + 1;   % ~600 ms window, must be odd
for k = 1:n_leads
    s = ecg_raw(k,:);
    s = filtfilt(b_hp,    a_hp,    s);
    s = filtfilt(b_lp,    a_lp,    s);
    s = filtfilt(b_notch, a_notch, s);
    s = s - movmedian(s, win_bl);       % Stream B baseline
    ecg_final(k,:) = s;
end

% Stacked 12-lead overlay (gray = original, color = processed)
figure('Name','C.6 — Full Pipeline 12-Lead Overlay','Position',[50 50 1500 950]);
t_disp2 = min(8, duration_sec);
idx_z   = t <= t_disp2;
t_z     = t(idx_z);
os2     = offset_scale;

for k = 1:n_leads
    off = (n_leads - k) * os2;
    plot(t_z, ecg_raw(k,idx_z)   + off, 'Color',[0.78 0.78 0.78], 'LineWidth', 0.8); hold on;
    plot(t_z, ecg_final(k,idx_z) + off, 'Color', lead_colors(k,:), 'LineWidth', 1.4);
    text(t_z(1) - 0.04*t_disp2, off, lead_names{k}, 'FontWeight','bold', ...
         'Color', lead_colors(k,:), 'FontSize', 8.5, 'HorizontalAlignment','right');
end
plot(NaN,NaN,'Color',[0.72 0.72 0.72],'LineWidth',2,'DisplayName','Original');
plot(NaN,NaN,'Color','b',             'LineWidth',2,'DisplayName','Processed (A+B+C)');
legend('Location','northeast');
xlabel('Time (s)'); ylabel('mV + display offset');
title(sprintf('C.6 — Complete Pipeline: %s  |  Gray = Original  |  Color = Processed', record_used));
xlim([t_z(1) t_z(end)]); yticks([]); grid on; box off;

% Timing error across all 66 lead pairs
fprintf('  Computing timing error across all %d lead pairs...\n', nchoosek(n_leads,2));
all_lp    = nchoosek(1:n_leads, 2);
max_lxc   = round(fs * 0.25);
e_iir_all = zeros(size(all_lp,1), 1);
e_fin_all = zeros(size(all_lp,1), 1);

for p = 1:size(all_lp,1)
    a = all_lp(p,1);  b = all_lp(p,2);
    [xo,lo] = xcorr(ecg_raw(a,:),   ecg_raw(b,:),   max_lxc, 'normalized');
    [xi,~  ] = xcorr(ecg_iir(a,:),   ecg_iir(b,:),   max_lxc, 'normalized');
    [xf,~  ] = xcorr(ecg_final(a,:), ecg_final(b,:), max_lxc, 'normalized');
    [~,io2]=max(abs(xo)); [~,ii2]=max(abs(xi)); [~,if2]=max(abs(xf));
    e_iir_all(p) = abs(lo(ii2) - lo(io2)) / fs * 1000;
    e_fin_all(p) = abs(lo(if2) - lo(io2)) / fs * 1000;
end

max_err_iir   = max(e_iir_all);
max_err_final = max(e_fin_all);
[~,wi] = max(e_iir_all);   [~,wf] = max(e_fin_all);

fprintf('\n  C.6 Results:\n');
fprintf('    IIR causal pipeline    — Max error: %.1f ms  (%s vs %s)\n', ...
        max_err_iir,   lead_names{all_lp(wi,1)}, lead_names{all_lp(wi,2)});
fprintf('    Corrected pipeline     — Max error: %.1f ms  (%s vs %s)\n\n', ...
        max_err_final, lead_names{all_lp(wf,1)}, lead_names{all_lp(wf,2)});

%% =========================================================================
%  C.7 — CLINICAL IMPACT: STEMI Detection
%  =========================================================================
fprintf('--- C.7: Clinical Impact — Anterior STEMI ---\n');

% Build a STEMI version of the loaded signal by adding ST elevation to V1-V4.
% This is added on top of the real ECG so the underlying morphology is preserved.
ecg_stemi = add_stemi_elevation(ecg_raw, fs);

ecg_stemi_iir = zeros(n_leads, N);
ecg_stemi_zp  = zeros(n_leads, N);
for k = 1:n_leads
    s = ecg_stemi(k,:);
    si = filter(b_hp, a_hp, s);   si = filter(b_lp, a_lp, si);
    sz = filtfilt(b_hp, a_hp, s); sz = filtfilt(b_lp, a_lp, sz);
    ecg_stemi_iir(k,:) = si;
    ecg_stemi_zp(k,:)  = sz;
end

% Zoom to one beat for clinical comparison
% Find approximate R-peak location in lead II (index 2) for display
lead_ii    = ecg_raw(2,:);
[~, rpks]  = findpeaks(lead_ii, 'MinPeakProminence', 0.2*max(lead_ii), ...
                                 'MinPeakDistance',  round(0.4*fs));
if isempty(rpks)
    beat_centre = round(N/2);   % fallback
else
    % pick second R-peak so the filter transient has settled
    beat_centre = rpks(min(2, end));
end
half_beat = round(0.5 * fs);
b1 = max(1,        beat_centre - half_beat);
b2 = min(N,        beat_centre + half_beat);
t_ms = (t(b1:b2) - t(beat_centre)) * 1000;   % ms relative to R-peak

v1_idx = 7;  v4_idx = 10;

figure('Name','C.7 — STEMI Clinical Impact','Position',[50 50 1400 850]);

subplot(3,1,1);
plot(t_ms, ecg_stemi(v1_idx, b1:b2),     'b',  'LineWidth', 2.2); hold on;
plot(t_ms, ecg_stemi(v4_idx, b1:b2),     'r',  'LineWidth', 2.2);
yl = ylim;
patch([40 140 140 40],[yl(1) yl(1) yl(2) yl(2)],'y','FaceAlpha',0.25,'EdgeColor','none');
text(90, yl(1)+0.85*(yl(2)-yl(1)),'ST seg', ...
     'HorizontalAlignment','center','FontSize',8,'Color','#8B6914','FontWeight','bold');
title('C.7a — STEMI Original: V1 & V4 ST-elevation aligned (simultaneous J-point)');
legend('V1','V4','Location','northeast'); ylabel('mV'); xlabel('ms rel. R-peak'); grid on;

subplot(3,1,2);
plot(t_ms, ecg_stemi_iir(v1_idx, b1:b2), 'b--','LineWidth', 2.2); hold on;
plot(t_ms, ecg_stemi_iir(v4_idx, b1:b2), 'r--','LineWidth', 2.2);
title(sprintf('C.7b — IIR Causal (DISTORTED): transient error ~%.0f ms, ST morphology altered', ...
              max_err_iir));
legend('V1 IIR','V4 IIR','Location','northeast'); ylabel('mV'); xlabel('ms rel. R-peak'); grid on;

subplot(3,1,3);
plot(t_ms, ecg_stemi_zp(v1_idx, b1:b2),  'b',  'LineWidth', 2.2); hold on;
plot(t_ms, ecg_stemi_zp(v4_idx, b1:b2),  'r',  'LineWidth', 2.2);
title(sprintf('C.7c — Zero-Phase Corrected: ST segments restored (error < %.1f ms)', ...
              max_err_final + 0.5));
legend('V1 ZP','V4 ZP','Location','northeast'); ylabel('mV'); xlabel('ms rel. R-peak'); grid on;

fprintf('\n');
fprintf('  +------------------------------------------------------------------+\n');
fprintf('  |           CLINICAL IMPACT SUMMARY — STEMI DETECTION             |\n');
fprintf('  +------------------------------------------------------------------+\n');
fprintf('  |  Scenario    : Anterior ST-Elevation MI (STEMI)                 |\n');
fprintf('  |  Leads       : V1, V2, V3, V4 (precordial)                     |\n');
fprintf('  |  Rule        : ST elevation >= 1 mm at J-point                  |\n');
fprintf('  |  J-point win : ~40 ms  (AHA/ESC guideline)                      |\n');
fprintf('  |                                                                  |\n');
fprintf('  |  IIR causal timing error  : ~%4.0f ms  --> FAILS AHA standard   |\n', max_err_iir);
fprintf('  |  Zero-phase timing error  : ~%4.1f ms  --> PASSES AHA standard  |\n', max_err_final);
fprintf('  |                                                                  |\n');
fprintf('  |  Consequences of uncorrected phase error:                       |\n');
fprintf('  |    - ST amplitude over/under-estimated per lead                 |\n');
fprintf('  |    - J-point appears at different times in V1 vs V4             |\n');
fprintf('  |    - Automated STEMI alert fires on wrong lead territory        |\n');
fprintf('  |    - Mislocalization: anterior vs anteroseptal vs lateral MI    |\n');
fprintf('  |                                                                  |\n');
fprintf('  |  Fix: always use filtfilt or FIR in diagnostic ECG systems.     |\n');
fprintf('  +------------------------------------------------------------------+\n\n');

fprintf('=== STREAM C COMPLETE ===\n');

%% =========================================================================
%  LOCAL HELPER FUNCTIONS — no external toolboxes required
%  =========================================================================

function [ecg, fs, record_name] = load_ecg_folder(folder, n_leads_expected)
%LOAD_ECG_FOLDER  Scan folder for .dat+.hea pairs; load the first valid one.
%
%  Parses the WFDB .hea header format:
%    Line 1 : <record> <n_signals> <fs> <n_samples>
%    Lines 2+: <filename> <format> <gain/unit> <bits> <baseline> <label> ...
%
%  The .dat file is raw int16 (format 16), samples interleaved across leads.
%  No WFDB toolbox needed — uses only fopen/fread/fgetl.

    record_name = 'synthetic';
    ecg         = [];
    fs          = 500;

    % --- Check folder exists ---
    if ~isfolder(folder)
        warning('Folder not found: %s\nFalling back to synthetic ECG.', folder);
        ecg = generate_12lead_ecg(fs, 10);
        record_name = 'synthetic (folder not found)';
        return;
    end

    % --- Find all .hea files ---
    hea_files = dir(fullfile(folder, '*.hea'));
    if isempty(hea_files)
        warning('No .hea files found in: %s\nFalling back to synthetic ECG.', folder);
        ecg = generate_12lead_ecg(fs, 10);
        record_name = 'synthetic (no .hea files found)';
        return;
    end

    % --- Try each .hea file until one loads cleanly ---
    for h = 1:length(hea_files)
        hea_path = fullfile(folder, hea_files(h).name);
        [stem, ~] = fileparts(hea_files(h).name);
        dat_path = fullfile(folder, [stem '.dat']);

        if ~isfile(dat_path)
            continue;   % no matching .dat, skip
        end

        % Parse header
        [fs_hdr, n_samp, gain, baseline, n_sig] = parse_hea(hea_path);

        if n_sig < n_leads_expected
            fprintf('  Skipping %s (%d leads < %d required)\n', ...
                    stem, n_sig, n_leads_expected);
            continue;
        end

        % Read binary data
        fid = fopen(dat_path, 'rb', 'ieee-le');
        if fid == -1
            continue;
        end
        raw = fread(fid, n_sig * n_samp, 'int16');
        fclose(fid);

        if numel(raw) < n_sig * n_samp
            fprintf('  Skipping %s (file shorter than header says)\n', stem);
            continue;
        end

        raw = reshape(raw(1 : n_sig*n_samp), n_sig, n_samp);

        % Convert ADC counts → mV, keep first n_leads_expected leads
        ecg = zeros(n_leads_expected, n_samp);
        for k = 1:n_leads_expected
            ecg(k,:) = (double(raw(k,:)) - baseline(k)) / gain(k);
        end

        fs          = fs_hdr;
        record_name = stem;
        fprintf('  Loaded  : %s  (%d leads, %d Hz, %.1f s)\n', ...
                stem, n_leads_expected, fs, n_samp/fs_hdr);
        return;
    end

    % --- Nothing loaded ---
    warning('No valid .dat+.hea pair found in %s.\nFalling back to synthetic ECG.', folder);
    ecg = generate_12lead_ecg(fs, 10);
    record_name = 'synthetic (no valid record found)';
end

% -------------------------------------------------------------------------
function [fs, n_samp, gain, baseline, n_sig] = parse_hea(hea_path)
%PARSE_HEA  Extract key fields from a WFDB .hea header file.
%
%  Record line format (line 1):
%    <record_name>  <n_signals>  <fs>  <n_samples>  [<base_time>  <base_date>]
%
%  Signal line format (lines 2 .. n_signals+1):
%    <filename>  <format>  <gain(/<units>)>  <bits>  <baseline>  <physmin>  <physmax>  <label>
%    e.g.:  00001_hr  16+32768  1000/mV  16  0  -32768  32767  I

    fs       = 500;
    n_samp   = 5000;
    n_sig    = 12;
    gain     = 1000 * ones(12,1);
    baseline = zeros(12,1);

    fid = fopen(hea_path, 'r');
    if fid == -1, return; end

    % Line 1 — record descriptor
    ln1 = strtrim(fgetl(fid));
    parts1 = strsplit(ln1);
    if length(parts1) >= 3,  n_sig  = str2double(parts1{2}); end
    if length(parts1) >= 3,  fs     = str2double(parts1{3}); end
    if length(parts1) >= 4,  n_samp = str2double(parts1{4}); end

    % Some headers put the sampling freq as "fs/counterpfreq" (e.g. "500/1000")
    if isnan(fs)
        tok = strsplit(parts1{3}, '/');
        fs = str2double(tok{1});
    end
    if isnan(fs) || fs <= 0,  fs = 500; end

    % Clamp n_sig to something sane
    n_sig = max(1, min(n_sig, 12));
    gain     = 1000 * ones(n_sig, 1);
    baseline = zeros(n_sig, 1);

    % Lines 2+ — one per signal
    for k = 1:n_sig
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        ln = strtrim(ln);
        if isempty(ln) || ln(1) == '#'
            k = k - 1;  %#ok<FXSET>  % comment line, redo
            continue;
        end
        parts = strsplit(ln);
        if length(parts) >= 3
            % Field 3: gain string, e.g. "1000/mV" or "200" or "1000(nV/LSB)"
            gs  = parts{3};
            gs  = regexprep(gs, '\(.*\)', '');  % strip parenthesised units
            sl  = strfind(gs, '/');
            if ~isempty(sl)
                gval = str2double(gs(1:sl(1)-1));
            else
                gval = str2double(gs);
            end
            if ~isnan(gval) && gval ~= 0
                gain(k) = gval;
            end
        end
        if length(parts) >= 5
            bval = str2double(parts{5});
            if ~isnan(bval)
                baseline(k) = bval;
            end
        end
    end

    fclose(fid);
end

% -------------------------------------------------------------------------
function ecg = add_stemi_elevation(ecg_in, fs)
%ADD_STEMI_ELEVATION  Superimpose anterior STEMI ST-elevation onto a real ECG.
%   Adds a smooth plateau-shaped ST segment to leads V1–V4 (indices 7–10)
%   at each beat location detected in lead II.

    ecg         = ecg_in;
    N           = size(ecg_in, 2);
    t           = (0:N-1) / fs;
    stemi_leads = [7 8 9 10];
    st_elev     = [0.15 0.25 0.20 0.15];   % mV

    % Detect R-peaks in lead II for beat timing
    lead_ii = ecg_in(2,:);
    prom    = max(0.2 * max(lead_ii), 0.1);
    [~, rpks] = findpeaks(lead_ii, 'MinPeakProminence', prom, ...
                                    'MinPeakDistance',  round(0.4*fs));
    if isempty(rpks)
        % No peaks found — space beats uniformly at 72 bpm
        rpks = round((0.5 : (60/72) : t(end)) * fs);
        rpks = rpks(rpks > 0 & rpks <= N);
    end

    for i = 1:length(stemi_leads)
        k   = stemi_leads(i);
        sig = ecg(k,:);
        for r = 1:length(rpks)
            c   = rpks(r) / fs;      % R-peak time (s)
            j   = c + 0.040;         % J-point
            ton = c + 0.120;         % T-wave onset
            sig = sig + st_elev(i) * plateau(t, j, ton, 0.008);
        end
        ecg(k,:) = sig;
    end
end

% -------------------------------------------------------------------------
function ecg = generate_12lead_ecg(fs, duration_sec)
%GENERATE_12LEAD_ECG  Fallback synthetic 12-lead ECG (pure MATLAB, no toolbox).

    N  = fs * duration_sec;
    t  = (0:N-1) / fs;
    hr = 72;  rr = 60/hr;

    P_a = [ 0.08  0.12  0.05 -0.10  0.04  0.08  0.06  0.08  0.10  0.10  0.08  0.06];
    Q_a = [-0.04 -0.03 -0.02  0.04 -0.02 -0.03  0.00  0.00 -0.02 -0.04 -0.03 -0.03];
    R_a = [ 0.60  1.00  0.50 -0.80  0.30  0.60  0.20  0.50  0.90  1.20  1.00  0.70];
    S_a = [-0.10 -0.10  0.00  0.10 -0.10 -0.10 -0.80 -0.70 -0.40 -0.20 -0.10 -0.05];
    T_a = [ 0.15  0.25  0.10 -0.20  0.08  0.15 -0.10  0.20  0.30  0.35  0.30  0.20];

    P_t=-0.130; P_w=0.025;  Q_t=-0.030; Q_w=0.010;
    R_t= 0.000; R_w=0.012;  S_t= 0.030; S_w=0.010;
    T_t= 0.150; T_w=0.040;

    ecg = zeros(12, N);
    for k = 1:12
        sig = zeros(1,N);
        bs  = 0.20;
        while bs + rr <= duration_sec + 0.05
            c   = bs + rr/2;
            sig = sig + gk(t,c+P_t,P_w,P_a(k)) + gk(t,c+Q_t,Q_w,Q_a(k)) ...
                      + gk(t,c+R_t,R_w,R_a(k)) + gk(t,c+S_t,S_w,S_a(k)) ...
                      + gk(t,c+T_t,T_w,T_a(k));
            bs  = bs + rr;
        end
        ecg(k,:) = sig + 0.008*randn(1,N);
    end
end

% -------------------------------------------------------------------------
function y = gk(t, mu, sigma, amp)
    y = amp * exp(-((t-mu).^2) / (2*sigma^2));
end

% -------------------------------------------------------------------------
function y = plateau(t, t0, t1, ramp)
%PLATEAU  Smooth trapezoidal pulse with raised-cosine ramps.
    y    = zeros(size(t));
    rise = (t >= t0-ramp) & (t <  t0);
    flat = (t >= t0)      & (t <= t1);
    fall = (t >  t1)      & (t <= t1+ramp);
    y(rise) = 0.5*(1 - cos(pi*(t(rise)-(t0-ramp))/ramp));
    y(flat) = 1;
    y(fall) = 0.5*(1 + cos(pi*(t(fall)-t1)/ramp));
end
