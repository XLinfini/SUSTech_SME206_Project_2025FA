clear; close all; clc;

function detection_mask = cfar_detector_2d(rdm, train_cells_range, guard_cells_range, train_cells_doppler, guard_cells_doppler, snr_offset_db)

[num_range_cells, num_doppler_cells] = size(rdm);
detection_mask = zeros(num_range_cells, num_doppler_cells);

snr_offset = 10^(snr_offset_db / 10);

win_size_range = 2 * (train_cells_range + guard_cells_range) + 1;
win_size_doppler = 2 * (train_cells_doppler + guard_cells_doppler) + 1;

rdm_padded = padarray(rdm, [train_cells_range + guard_cells_range, train_cells_doppler + guard_cells_doppler], 'replicate');

for r = 1:num_range_cells
    for d = 1:num_doppler_cells
        
        cut_r = r + train_cells_range + guard_cells_range;
        cut_d = d + train_cells_doppler + guard_cells_doppler;
        
        r_start = cut_r - (train_cells_range + guard_cells_range);
        r_end = cut_r + (train_cells_range + guard_cells_range);
        d_start = cut_d - (train_cells_doppler + guard_cells_doppler);
        d_end = cut_d + (train_cells_doppler + guard_cells_doppler);
        training_window = rdm_padded(r_start:r_end, d_start:d_end);
        
        g_r_start = train_cells_range + 1;
        g_r_end = g_r_start + 2 * guard_cells_range;
        g_d_start = train_cells_doppler + 1;
        g_d_end = g_d_start + 2 * guard_cells_doppler;
        training_window(g_r_start:g_r_end, g_d_start:g_d_end) = 0;
        
        num_train_cells = numel(training_window) - (2*guard_cells_range+1)*(2*guard_cells_doppler+1);
        
        noise_level = sum(training_window(:)) / num_train_cells;
        
        threshold = noise_level * snr_offset;
        
        if rdm(r, d) > threshold
            detection_mask(r, d) = 1;
        end
    end
end

end

fc = 77e9;
B = 150e6;
Tchirp = 8e-6;
Nr = 1024;
Nd = 128;
c = 3e8;

K = B / Tchirp;
lambda = c / fc;
fs = Nr / Tchirp;

targets = struct();
targets(1).R0 = 80;   targets(1).v0 = 10;
targets(2).R0 = 120;  targets(2).v0 = -20;
targets(3).R0 = 180;  targets(3).v0 = 30;

t_fast = (0:Nr-1) / fs;
t_slow = (0:Nd-1) * Tchirp;

if_signal = zeros(Nr, Nd);

for m = 1:length(targets)
    for i = 1:Nd
        t_current_slow = t_slow(i);
        current_R = targets(m).R0 + targets(m).v0 * t_current_slow;
        tau = 2 * current_R / c;
        
        phase = 2*pi * ( (K*tau) * t_fast.' + (2*fc*targets(m).v0/c) * t_current_slow );
        if_signal(:, i) = if_signal(:, i) + cos(phase);
    end
end
noise_power = 0.05;
if_signal = if_signal + sqrt(noise_power/2) * (randn(size(if_signal)) + 1i*randn(size(if_signal)));


range_fft = fft(if_signal, Nr);
doppler_fft = fftshift(fft(range_fft.', Nd), 1).';
rdm = abs(doppler_fft);
rdm = rdm(1:Nr/2, :);

train_cells_range = 8;
guard_cells_range = 4;
train_cells_doppler = 4;
guard_cells_doppler = 2;
snr_offset_db = 15;

detection_mask = cfar_detector_2d(rdm, train_cells_range, guard_cells_range, train_cells_doppler, guard_cells_doppler, snr_offset_db);

figure;
subplot(1, 2, 1);
range_resolution = c / (2 * B);
velocity_resolution = lambda / (2 * Nd * Tchirp);
range_axis = (0:Nr/2-1) * range_resolution;
velocity_axis = (-Nd/2 : Nd/2-1) * velocity_resolution;
imagesc(velocity_axis, range_axis, 10*log10(rdm));
colorbar;
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map');
axis xy;

subplot(1, 2, 2);
imagesc(velocity_axis, range_axis, detection_mask);
colorbar;
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('CFAR Detection Mask');
axis xy;

[cand_range_idx, cand_doppler_idx] = find(detection_mask);

detected_range_idx = [];
detected_doppler_idx = [];

for i = 1:length(cand_range_idx)
    r = cand_range_idx(i);
    d = cand_doppler_idx(i);
    
    r_start = max(1, r-1);
    r_end = min(size(rdm, 1), r+1);
    d_start = max(1, d-1);
    d_end = min(size(rdm, 2), d+1);
    
    neighborhood = rdm(r_start:r_end, d_start:d_end);
    
    if rdm(r, d) >= max(neighborhood(:))
        detected_range_idx = [detected_range_idx; r];
        detected_doppler_idx = [detected_doppler_idx; d];
    end
end

peak_indices = unique([detected_range_idx, detected_doppler_idx], 'rows');
detected_range_idx = peak_indices(:, 1);
detected_doppler_idx = peak_indices(:, 2);

fprintf('--- Multi-Target Detection Results ---\n');
fprintf('Ground Truth:\n');
for m = 1:length(targets)
    fprintf('  Target %d: Range=%.1f m, Velocity=%.1f m/s\n', m, targets(m).R0, targets(m).v0);
end

fprintf('Detections:\n');
if isempty(detected_range_idx)
    fprintf('  No targets detected.\n');
else
    for i = 1:length(detected_range_idx)
        measured_range = (detected_range_idx(i) - 1) * range_resolution;
        measured_velocity = (detected_doppler_idx(i) - (Nd/2 + 1)) * velocity_resolution;
        fprintf('  Detected Target %d: Range=%.2f m, Velocity=%.2f m/s\n', i, measured_range, measured_velocity);
    end
end
