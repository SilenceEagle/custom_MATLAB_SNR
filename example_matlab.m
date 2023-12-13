%% example usage of custom_matlab_snr function
clear
close all

rng(111); % reproducability

fs = 4e4;  % sampling rate in Hz
t = 0:1/fs:0.05;  % sampled time points

%% mono frequency sinusidal signal
disp('generating mono frequency signal...')
target_frequency = 200; % Hz
s1_clean = sin(2*pi*target_frequency*t);

% add noise
target_snr = 10;
s1 = awgn(s1_clean, target_snr, 'measured');

% display signal in time domain
figure;
plot(t, s1, 'k', 'DisplayName', 'noisy signal');
hold on
grid on
plot(t, s1_clean, 'r', 'LineWidth', 2, 'DisplayName', 'clean signal');
legend('Location', 'best');
xlabel('Time (s)')
ylabel('Value')

% count snr

% 1 harmonic, given target frequency
number_harmonic = 1;
% snr(s1, fs, number_harmonic);
custom_matlab_snr(s1, fs, number_harmonic, target_frequency);
[snr_s1, noise_power_s1] = custom_matlab_snr(s1, fs, number_harmonic, ...
   target_frequency);
fprintf('SNR with 1 harmonic: %.4f\n', snr_s1)

% 1 harmonic, given target freuency, exclude psd below 300 Hz when counting
% noise power
max_frequency_dc = target_frequency;  % Hz
custom_matlab_snr(s1, fs, number_harmonic, target_frequency, max_frequency_dc);
[snr_s1_dc, noise_power_s1_dc] = custom_matlab_snr(s1, fs, number_harmonic, ...
    target_frequency, max_frequency_dc);
fprintf('SNR with 1 harmonic and force DC: %.4f\n', snr_s1_dc)

% matlab snr from time domain results
snr_s1_t = snr(s1, s1-s1_clean);
fprintf('SNR from time domian: %.4f\n', snr_s1_t)

%% multi frequency
disp('generating dual frequency signal...')
target_frequency2 = 500;  % Hz
s2_clean = s1_clean + sin(2*pi*target_frequency2*t);

% add nosie
s2 = awgn(s2_clean, target_snr, 'measured');

% display signal in time domain
figure;
plot(t, s2, 'k', 'DisplayName', 'noisy signal');
hold on
grid on
plot(t, s2_clean, 'r', 'LineWidth', 2, 'DisplayName', 'clean signal');
legend('Location', 'best');
xlabel('Time (s)')
ylabel('Value')


% 1 harmonic, given target frequencies
disp('1 harmonic:')

% frequency 1
% snr(s2, fs, number_harmonic);
custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency, target_frequency2]); % signal power of 2nd frequency
    % will be excluded from counting noise power (By marking it as 2nd Harmonic)
[snr_s21, noise_power_s21] = custom_matlab_snr(s2, fs, number_harmonic,...
    [target_frequency, target_frequency2]);

% frequency 2
custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency2, target_frequency]);  % move target frequency to be first
[snr_s22, noise_power_s22] = custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency2, target_frequency]);

% overall snr
snr_s2 = mean([snr_s21, snr_s22]);
fprintf('SNR:\n\tOverall(mean): %.4f dB\n\t%d Hz: %.4f dB\n\t%d Hz: %.4f dB\n',...
    snr_s2, target_frequency, snr_s21, target_frequency2, snr_s22);
fprintf('\tOverall(sqrt sum square): %.4f dB\n',...
    sqrt((snr_s21^2+snr_s22^2)))

% 1 harmonic, given target freuency, exclude psd below 200 Hz when counting
% noise power
disp('1 harmonic, exclude psd below 200 Hz when counting noise power:')
max_frequency_dc = target_frequency;  % Hz

% frequency 1
custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency, target_frequency2], max_frequency_dc);
[snr_s23, noise_power_s23] = custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency, target_frequency2], max_frequency_dc);

% frequency 2
custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency2, target_frequency], max_frequency_dc);
[snr_s24, noise_power_s24] = custom_matlab_snr(s2, fs, number_harmonic, ...
    [target_frequency2, target_frequency], max_frequency_dc);

snr_s2_c = mean([snr_s23, snr_s24]);
fprintf('SNR:\n\tOverall(mean): %.4f dB\n\t%d Hz: %.4f dB\n\t%d Hz: %.4f dB\n',...
    snr_s2_c, target_frequency, snr_s23, target_frequency2, snr_s24);
fprintf('\tOverall(sqrt sum square): %.4f dB\n',...
    sqrt((snr_s23^2+snr_s24^2)))


% matlab snr from time domain results
snr_s2_t = snr(s2, s2-s2_clean);
fprintf('SNR from time domian: %.4f\n', snr_s2_t)