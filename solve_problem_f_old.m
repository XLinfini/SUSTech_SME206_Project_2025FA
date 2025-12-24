% SME206 Project: Problem f) - Real Data Analysis with Triangular Modulation

% Clear workspace, close figures, and clear command window
clear; close all; clc;

% --- 1. System Parameters ---
fc = 2.253e9;       % Carrier frequency (Hz)
B = 178e6;          % Bandwidth (Hz)
Tt = 1/65;          % Triangle wave period (s)
fs = 50e3;          % Sampling frequency (Hz)
c = 3e8;            % Speed of light (m/s)

% Calculated parameters
K = 2 * B / Tt;     % Triangle wave slope (Hz/s)
lambda = c / fc;

% --- 2. Load and Separate Data ---
fprintf('Loading real data file...\n');
data = readmatrix('fmcw_bonus_f(1).csv');
if_signal_raw = data(:, 1);
sync_signal = data(:, 2);
fprintf('Data loaded. Total samples: %d\n', length(if_signal_raw));

% --- 3. Synchronize and Segment Data ---
% Find the start of each chirp using the sync signal.
% A simple threshold-based method can find the rising edges.
sync_threshold = (max(sync_signal) + min(sync_signal)) / 2;
rising_edges_idx = find(diff(sync_signal > sync_threshold) == 1);

% Each period consists of an up-chirp and a down-chirp
num_periods = floor(length(rising_edges_idx));
fprintf('Detected %d full triangular periods.\n', num_periods);

% Initialize arrays to store results
measured_distances = zeros(num_periods, 1);
measured_velocities = zeros(num_periods, 1);
time_axis = (0:num_periods-1) * Tt;

% --- 4. Process Each Triangular Period ---
fprintf('Processing each period...\n');
% Loop to num_periods - 1 to avoid index out of bounds on the last element
for i = 1:(num_periods - 1)
    % Get indices for the current up and down chirp
    up_chirp_start = rising_edges_idx(i);
    up_chirp_end = round((rising_edges_idx(i) + rising_edges_idx(i+1)) / 2);
    down_chirp_start = up_chirp_end;
    down_chirp_end = rising_edges_idx(i+1);
    
    % Extract signal segments
    up_chirp_signal = if_signal_raw(up_chirp_start:up_chirp_end);
    down_chirp_signal = if_signal_raw(down_chirp_start:down_chirp_end);
    
    % --- Pre-processing ---
    % 1. Remove DC offset
    up_chirp_signal = up_chirp_signal - mean(up_chirp_signal);
    down_chirp_signal = down_chirp_signal - mean(down_chirp_signal);
    
    % Perform FFT on both segments
    % Use a large FFT point number for better frequency resolution
    N_fft = 2^nextpow2(length(up_chirp_signal)) * 4;
    
    up_fft = fft(up_chirp_signal, N_fft);
    down_fft = fft(down_chirp_signal, N_fft);

    % 2. Smooth the spectrum to reduce noise effects
    up_fft = movmean(up_fft, 5);
    down_fft = movmean(down_fft, 5);
    
    % Find peak frequencies
    [~, up_idx] = max(abs(up_fft(101:N_fft/2)));
    up_idx = up_idx + 100;
    [~, down_idx] = max(abs(down_fft(101:N_fft/2)));
    down_idx = down_idx + 100;

    % Store complex values of fft at peak freq
    up_complex_val = up_fft(up_idx);
    down_complex_val = down_fft(down_idx);

    % Construct freq axis to determine distance
    f_axis = fs * (0:N_fft/2-1) / N_fft;
    f_up = f_axis(up_idx);
    f_down = f_axis(down_idx);
    
    % Calculate beat and doppler frequencies
    % --- Range Calculation (Robust) ---
    f_beat = (f_up + f_down) / 2;
    measured_distances(i) = (c * f_beat) / (2 * K);

    % --- Velocity Calculation ---
    if i > 1
        % Calculate phase difference
        phase_diff_up = angle(up_complex_val * conj(prev_up_complex_val));
        phase_diff_down = angle(down_complex_val * conj(prev_down_complex_val));
        phase_diff = (phase_diff_up - phase_diff_down) / 2;
        
        % Calculate Doppler frequency and velocity
        f_doppler = phase_diff / (2 * pi * Tt);
        measured_velocities(i) = (f_doppler * lambda) / 2;
    end

    % Store this round's chirp indicies for next use
    prev_up_complex_val = up_fft(up_idx);
    prev_down_complex_val = down_fft(down_idx);
end
fprintf('Processing complete.\n');

% --- 5. Analyze and Visualize Results ---
% Filter out potential outliers for cleaner plots
% A simple median filter can work well for this
measured_distances = medfilt1(measured_distances, 5);
measured_velocities = medfilt1(measured_velocities, 5);

% Plot distance over time
figure;
subplot(2, 1, 1);
plot(time_axis, measured_distances);
grid on;
title('Distance Change of the Tester Over Time');
xlabel('Time (s)');
ylabel('Distance (m)');

% Plot velocity over time
subplot(2, 1, 2);
plot(time_axis, measured_velocities);
grid on;
title('Velocity Change of the Tester Over Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Display velocity range
min_velocity = min(measured_velocities);
max_velocity = max(measured_velocities);
fprintf('\n--- Analysis Results ---\n');
fprintf('Velocity range of the tester: [%.2f m/s, %.2f m/s]\n', min_velocity, max_velocity);