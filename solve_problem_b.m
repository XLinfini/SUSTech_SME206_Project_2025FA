clear; close all; clc;

fc = 77e9;
B = 150e6;
Tchirp = 8e-6;
K = B / Tchirp;
Nr = 1024;
Nd = 128;
c = 3e8;
lambda = c / fc;

data = load('fmcw_basic_b(1).mat');

x2_vec = data.x2;
x2_matrix = reshape(x2_vec, Nr, Nd);

doppler_fft = fftshift(fft2(x2_matrix, Nr, Nd), 2);

rdm = abs(doppler_fft);
rdm = rdm(1:Nr/2, :);

[max_val, max_idx_1d] = max(rdm(:));
[range_idx, doppler_idx] = ind2sub(size(rdm), max_idx_1d);

fs = Nr / Tchirp;
range_resolution = (c * fs) / (2 * K * Nr);
target_range = (range_idx - 1) * range_resolution;

doppler_freq_resolution = 1 / (Tchirp * Nd);
doppler_shift_freq = (doppler_idx - (Nd/2 + 1)) * doppler_freq_resolution;
target_velocity = (doppler_shift_freq * lambda) / 2;

range_axis = (0:Nr/2-1) * range_resolution;
velocity_axis = (-Nd/2 : Nd/2-1) * (lambda / (2 * Tchirp * Nd));

figure;
imagesc(velocity_axis, range_axis, 10*log10(rdm));
colorbar;
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title('Range-Doppler Map');
axis xy;

fprintf('Peak detected in RDM.\n');
fprintf('Calculated target distance: %.2f meters\n', target_range);
fprintf('Calculated target velocity: %.2f m/s\n', target_velocity);
fprintf('(Positive velocity means moving away from the radar)\n');
