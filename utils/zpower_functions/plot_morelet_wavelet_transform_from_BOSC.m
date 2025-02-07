%% Parameters
Fsample = 1000;         % Sampling frequency in Hz
wavenum = 6;            % Wave number for the Morlet wavelet
F = 1:100;              % Frequencies from 1 to 100 Hz (step = 1 Hz)

%% Generate a synthetic EEG signal (e.g., 5 seconds duration)
duration = length(eegsignal)/Fsample;                                   % Duration in seconds
t_signal = 0:1/Fsample:(duration - 1/Fsample);    % Time vector for EEG signal
% Synthetic signal: sum of 10 Hz and 20 Hz sine waves with added noise
% eegsignal = sin(2*pi*10*t_signal) + 0.5*sin(2*pi*20*t_signal) + 0.1*randn(size(t_signal));

%% Compute Time-Frequency Representation using Morlet Wavelets
% Initialize the time-frequency matrix B
B = zeros(length(F), length(eegsignal));

for idx = 1:length(F)
    freq = F(idx);
    % Standard deviation of the Gaussian envelope for the wavelet
    st = 1 / (2*pi*(freq/wavenum));
    % Normalization factor A
    A = 1 / sqrt(st * sqrt(pi));
    % Time vector for the wavelet: from -3.6*st to 3.6*st
    t_wave = -3.6*st : 1/Fsample : 3.6*st;
    % Generate the Morlet wavelet:
    % m(t) = A * exp(-t.^2/(2*st^2)) .* exp(1i*2*pi*freq*t)
    m = A * exp(-t_wave.^2 / (2*st^2)) .* exp(1i * 2*pi*freq*t_wave);
    
    % Convolve the EEG signal with the wavelet using the 'same' option
    y = conv(eegsignal, m, 'same');
    
    % Compute the power (squared amplitude) and store it in B
    B(idx, :) = abs(y).^2;
end

%% Plotting
figure;

% Left Panel: Example Morlet Wavelet at 10 Hz
f_example = 10;    % Example frequency for illustration
st_example = 1 / (2*pi*(f_example/wavenum));
A_example = 1 / sqrt(st_example * sqrt(pi));
t_wave_example = -3.6*st_example : 1/Fsample : 3.6*st_example;
m_example = A_example * exp(-t_wave_example.^2 / (2*st_example^2)) .* exp(1i * 2*pi*f_example*t_wave_example);

subplot(1, 2, 1);
plot(t_wave_example, real(m_example), 'b-', 'LineWidth', 1.5); hold on;
plot(t_wave_example, imag(m_example), 'r-', 'LineWidth', 1.5);
plot(t_wave_example, abs(m_example), 'k--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Morlet Wavelet at 10 Hz');
legend('Real Part', 'Imaginary Part', 'Envelope');
grid on;

% Right Panel: Time-Frequency Representation (Spectrogram)
subplot(1, 2, 2);
imagesc(t_signal, F, B);  % X-axis: time, Y-axis: frequency
axis xy;                % Ensure that low frequencies are at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Time-Frequency Representation');
colorbar;
colormap jet;
