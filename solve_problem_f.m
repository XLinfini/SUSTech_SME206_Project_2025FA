clear; close all; clc;

fc = 2.253e9;
B = 178e6;
Tt = 1/65;
fs = 50e3;
c = 3e8;

T_up = Tt / 2;
K = B / T_up;
lambda = c / fc;

data = readmatrix('fmcw_bonus_f(1).csv');
if_signal_raw = data(:, 1);
sync_signal = data(:, 2);

thresh = (max(sync_signal) + min(sync_signal)) / 2;
binary_sync = sync_signal > thresh;
rising_edges = find(diff(binary_sync) == 1);

samples_per_up = round(fs * T_up);
valid_chirps = [];
chirp_time_axis = [];

for i = 1:length(rising_edges)-1
    idx_start = rising_edges(i);
    if idx_start + samples_per_up < length(if_signal_raw)
        seg = if_signal_raw(idx_start : idx_start + samples_per_up - 1);
        valid_chirps = [valid_chirps, seg];
        chirp_time_axis = [chirp_time_axis, (i-1)*Tt];
    end
end

[Nr, Nd] = size(valid_chirps);

mean_chirp = mean(valid_chirps, 2);
chirps_mti = valid_chirps - mean_chirp;

win = hamming(Nr);
chirps_win = chirps_mti .* win;
N_fft = 2048;
range_profile_complex = fft(chirps_win, N_fft, 1);
range_profile_abs = abs(range_profile_complex(1:N_fft/2, :));

f_axis = (0:N_fft/2-1) * (fs / N_fft);
range_axis = (c * f_axis) / (2 * K);

range_profile_abs = movmean(range_profile_abs, 5, 2);

detected_range = zeros(1, Nd);
detected_velocity = zeros(1, Nd);

search_min = 15; search_max = 25;
[~, idx_15m] = min(abs(range_axis - search_min));
[~, idx_25m] = min(abs(range_axis - search_max));

first_avg_spec = mean(range_profile_abs(:, 1:10), 2);
roi_spec = zeros(size(first_avg_spec));
roi_spec(idx_15m:idx_25m) = first_avg_spec(idx_15m:idx_25m);
[~, current_idx] = max(roi_spec);

window_m = 3.0;
win_bins = round(window_m / (range_axis(2) - range_axis(1)));

for i = 1:Nd
    idx_start = max(1, current_idx - win_bins);
    idx_end = min(N_fft/2, current_idx + win_bins);
    
    curr_spec = range_profile_abs(:, i);
    window_spec = curr_spec(idx_start:idx_end);
    
    [~, local_peak_idx] = max(window_spec);
    
    current_idx = idx_start + local_peak_idx - 1;
    detected_range(i) = range_axis(current_idx);
    
    if i > 1
        curr_phasor = range_profile_complex(current_idx, i);
        prev_phasor = range_profile_complex(current_idx, i-1);
        
        d_phi = angle(curr_phasor * conj(prev_phasor));
        v_inst = (lambda * d_phi) / (4 * pi * Tt);
        
        if abs(v_inst) < 3.0
            detected_velocity(i) = v_inst;
        else
            detected_velocity(i) = detected_velocity(i-1);
        end
    end
end

detected_range = medfilt1(detected_range, 10);
detected_velocity = medfilt1(detected_velocity, 15);

figure('Position', [100, 100, 1000, 600]);

subplot(2, 2, [1, 3]);
range_img = 20*log10(range_profile_abs);
valid_max = max(range_img(idx_15m:end, :), [], 'all');
imagesc(chirp_time_axis, range_axis, range_img);
axis xy;
colormap('jet');
colorbar;
caxis([valid_max-20, valid_max]);
hold on;
plot(chirp_time_axis, detected_range, 'w--', 'LineWidth', 2);
hold off;
title('Range-Time Map (Strict Window Tracking)');
xlabel('Time (s)');
ylabel('Range (m)');
ylim([0, 70]);

subplot(2, 2, 2);
plot(chirp_time_axis, detected_range, 'b', 'LineWidth', 1.5);
grid on;
title('Distance Trace');
xlabel('Time (s)');
ylabel('Distance (m)');
ylim([0, 70]);

subplot(2, 2, 4);
plot(chirp_time_axis, detected_velocity, 'r', 'LineWidth', 1.5);
grid on;
title('Velocity Trace');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
ylim([-2, 2]);

fprintf('--- Final Results ---\n');
fprintf('Tracked Distance: %.1f m to %.1f m\n', min(detected_range), max(detected_range));
