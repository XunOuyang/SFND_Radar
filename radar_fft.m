%   Define a signal: amplitude = A, frequency = f
%       signal = A*cos(2*pi*f*t);
%   Run the fft for dimension of samples N
%       signal_fft = fft(signal, N);
%   The output of the FFT processing of a signal is a complex number a+j*b
%       Magnitude:
%           signal_fft = abs(signal_fft);
%   

clear all;
clc;

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

% TODO: Form a signal containing a 77 Hz sinusoid of amplitude 0.7 and a 43Hz 
% sinusoid of amplitude 2.
A1 = 0.7
f1 = 77

S1 = A1*cos(2*pi*f1*t);

A2 = 2.0
f2 = 43

S2 = A2*cos(2*pi*f2*t);

S = S1 + S2;

% Corrupt the signal with noise 
X = S + 2*randn(size(t));

% Plot the noisy signal in the time domain. It is difficult to identify the frequency components by looking at the signal X(t). 
plot(1000*t(1:50) ,X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

% TODO : Compute the Fourier transform of the signal. 

signal_fft = fft(X);
signal_fft_Mag = abs(signal_fft);
N = length(signal_fft_Mag);
P1 = signal_fft_Mag;

% TODO : Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.


% Plotting
f = Fs*(0:(L/2))/L;
plot(f,P1(1:length(f))) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')