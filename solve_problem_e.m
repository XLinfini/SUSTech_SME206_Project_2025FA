clear; close all; clc;

fc = 77e9;
B = 150e6;
Tchirp = 8e-6;
Nr = 1024;
Nd = 128;
c = 3e8;

Nx = 8;
lambda = c / fc;
d = lambda / 2;

K = B / Tchirp;
fs = Nr / Tchirp;

targets = struct();
targets(1).R0 = 100;  targets(1).v0 = 15;   targets(1).theta_deg = -20;
targets(2).R0 = 150;  targets(2).v0 = -25;  targets(2).theta_deg = 30;

t_fast = (0:Nr-1) / fs;
t_slow = (0:Nd-1) * Tchirp;

radar_cube = zeros(Nr, Nd, Nx);

for m = 1:length(targets)
    for n = 1:Nx
        for i = 1:Nd
            t_current_slow = t_slow(i);
            current_R = targets(m).R0 + targets(m).v0 * t_current_slow;
            tau = 2 * current_R / c;
            
            theta_rad = deg2rad(targets(m).theta_deg);
            phase_shift_angle = 2*pi * d * (n-1) * sin(theta_rad) / lambda;
            
            phase = 2*pi * ( (K*tau) * t_fast.' + (2*fc*targets(m).v0/c) * t_current_slow );
            
            radar_cube(:, i, n) = radar_cube(:, i, n) + cos(phase + phase_shift_angle);
        end
    end
end
noise_power = 0.1;
radar_cube = radar_cube + sqrt(noise_power/2) * (randn(size(radar_cube)) + 1i*randn(size(radar_cube)));

range_fft_cube = fft(radar_cube, Nr, 1);

doppler_fft_cube = fftshift(fft(range_fft_cube, Nd, 2), 2);

rdm = squeeze(sum(abs(doppler_fft_cube), 3));
rdm = rdm(1:Nr/2, :);

detection_mask = cfar_detector_2d(rdm, 8, 4, 4, 2, 15);
peak_mask = imregionalmax(rdm .* detection_mask);
[detected_range_idx, detected_doppler_idx] = find(peak_mask);

fprintf('--- 3D Localization Results ---\n');
fprintf('Ground Truth:\n');
for m = 1:length(targets)
    fprintf('Target %d: R=%.1f m, v=%.1f m/s, theta=%.1f deg\n', m, targets(m).R0, targets(m).v0, targets(m).theta_deg);
end

fprintf('Detections:\n');
if isempty(detected_range_idx)
    fprintf('  No targets detected.\n');
else
    for i = 1:length(detected_range_idx)
        r_idx = detected_range_idx(i);
        d_idx = detected_doppler_idx(i);
        
        antenna_vec = squeeze(doppler_fft_cube(r_idx, d_idx, :));
        
        angle_fft = fftshift(fft(antenna_vec, Nx));
        
        [~, angle_idx] = max(abs(angle_fft));
        
        angle_sin = (angle_idx - (Nx/2 + 1)) * (2/Nx);
        measured_angle_deg = asind(angle_sin);
        
        range_resolution = c / (2 * B);
        velocity_resolution = lambda / (2 * Nd * Tchirp);
        measured_range = (r_idx - 1) * range_resolution;
        measured_velocity = (d_idx - (Nd/2 + 1)) * velocity_resolution;
        
        fprintf('Target %d: R=%.2f m, v=%.2f m/s, theta=%.2f deg\n', i, measured_range, measured_velocity, measured_angle_deg);
    end
end
