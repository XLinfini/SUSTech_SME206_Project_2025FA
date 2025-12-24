clear; close all; clc;

fc = 77e9;
B = 150e6;
Tchirp = 8e-6;
K = B / Tchirp;
Nr = 1024;
c = 3e8;

data = load('fmcw_basic_a(1).mat');
x1 = data.x1;

fs = Nr / Tchirp;

Y = fft(x1, Nr);

P2 = abs(Y / Nr);
P1 = P2(1:Nr/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs * (0:(Nr/2)) / Nr;

[peak_amplitude, idx] = max(P1);
fm = f(idx);

range = (c * fm) / (2 * K);

figure;
plot(f / 1e6, P1);
grid on;
title('Single-Sided Amplitude Spectrum of IF Signal');
xlabel('Frequency (MHz)');
ylabel('Amplitude');
hold on;
plot(fm / 1e6, peak_amplitude, 'r*', 'MarkerSize', 10);
legend('Spectrum', sprintf('Detected Peak (%.2f MHz)', fm/1e6));
hold off;

fprintf('Detected IF frequency (fm): %.2f MHz\n', fm / 1e6);
fprintf('Calculated target distance: %.2f meters\n', range);
