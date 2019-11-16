% Implement 1D CFAR using lagging cells on the given noise and target scenario.

% Close and delete all currently open figures
close all;
clear all;

% Data_points
Ns = 1000;

% Generate random noise
s=randn(Ns,1);

%Targets location. Assigning bin 100, 200, 300 and 700 as Targets with the amplitudes of 8, 9, 4, 11.
s([100 ,200, 300, 700])=[8 9 4 11];

%plot the output
plot(s);

% TODO: Apply CFAR to detect the targets by filtering the noise.

% 1. Define the following:
% 1a. Training Cells
% 1b. Guard Cells 
T = 3;
G = 1;

% Offset : Adding room above noise threshold for desired SNR 
% CA-CFAR3.5 to 4 may miss target or see the noise as target
offset=2.8;

% Vector to hold threshold values 
%threshold_cfar = zeros(1, Ns);
threshold_cfar = [];

%Vector to hold final signal after thresholding
signal_cfar = [];

% 2. Slide window across the signal length
windows_size = 0;
total_sum = 0;
avg = 0;
alpha = 0.9;
    

for i = (G+T+1):(Ns-(G+T))     
    %% 2. - 5. Determine the noise threshold by measuring it within the training cells
    cell_left = s(i-G-T:i-G-1);
    cell_right = s(i+G+1:i+G+T);
    % CA-CFAR
    Z = alpha * (sum(cell_left)+sum(cell_right))./(2*T);
    % GO-CFAR
    Z = max(max(cell_left), max(cell_right));
        
    threshold_cfar = [threshold_cfar, {Z}];
    % 6. Measuring the signal within the CUT
    
    % 8. Filter the signal above the threshold
    if s(i) - offset >= Z
        signal = s(i);
    else
        signal = 0;
    end
    signal_cfar = [signal_cfar, {signal}];
end


% plot the filtered signal
plot (cell2mat(signal_cfar),'g--');

% plot original sig, threshold and filtered signal within the same figure.
figure,plot(s);
hold on,plot(cell2mat(circshift(threshold_cfar,G)),'r--','LineWidth',2)
hold on, plot (cell2mat(circshift(signal_cfar,(T+G))),'g--','LineWidth',4);
legend('Signal','CFAR Threshold','detection')