%%% Scatter Plot and Histogram @ SNR = 3dB and SNR = 10dB
%%% Written by P.-J. Su 2022/8/30

clc;    %% Clear command window
clear;  %% Remove items from workspace

Eb = 1;     %% Bit energy (J)
T = 1;      %% Bit interval (s)
num_bits = 10^6;    %% Number of data bit
    
%%%  AWGN channel (SNR = 3dB) %%%
SNR_3dB = 10.^(3/10);
N0_3dB = 2*Eb/SNR_3dB;  %% SNR = 2*Eb/N0
sigma_3dB = sqrt(N0_3dB*T);

% When data bit d = +1 is transmitted (Source bit: 1)
z_1_3dB = sqrt(2*Eb*T) + randn(1, num_bits)*sigma_3dB;
% When data bit d = -1 is transmitted (Source bit: 0)
z_0_3dB = -sqrt(2*Eb*T) + randn(1, num_bits)*sigma_3dB;

%%%  AWGN channel (SNR = 10dB) %%%
SNR_10dB = 10.^(10/10);
N0_10dB = 2*Eb/SNR_10dB;    %% SNR = 2*Eb/N0
sigma_10dB = sqrt(N0_10dB*T);

% When data bit d = +1 is transmitted (Source bit: 1)
z_1_10dB = sqrt(2*Eb*T) + randn(1, num_bits)*sigma_10dB;
% When data bit d = -1 is transmitted (Source bit: 0)
z_0_10dB = -sqrt(2*Eb*T) + randn(1, num_bits)*sigma_10dB;

%%% Plot the scatter plot and histogram @ SNR = 3dB %%%
clf;
y_Idx = zeros(1, num_bits);
scatter(z_1_3dB,y_Idx,'filled');
hold on;
scatter(z_0_3dB,y_Idx,'filled','d');
title('Scatter Plot @ SNR = 3dB');
yticks([]);
legend('d = +1','d = -1');
xlabel('z');

figure;
nbins = 50;
h1 = histogram(z_1_3dB,nbins,'Normalization','pdf');
hold on;
h2 = histogram(z_0_3dB,nbins,'Normalization','pdf');
ylim([0 1]);
title('Histogram @ SNR = 3dB');
legend('\it f\rm_z(z|d = +1)','\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('\it f\rm(z)');

%%% Plot the scatter plot and histogram @ SNR = 10dB %%%
figure;
scatter(z_1_10dB,y_Idx,'filled');
hold on;
scatter(z_0_10dB,y_Idx,'filled','d');
title('Scatter Plot @ SNR = 10dB');
yticks([]);
legend('d = +1','d = -1');
xlabel('z');

figure;
h3 = histogram(z_1_10dB,nbins,'Normalization','pdf');
hold on;
h4 = histogram(z_0_10dB,nbins,'Normalization','pdf');
x_axis_max = sqrt(2*Eb*T) + 6*sigma_10dB;
x_axis_min = -sqrt(2*Eb*T) - 6*sigma_10dB;
axis([x_axis_min x_axis_max 0 1]);
title('Histogram @ SNR = 10dB');
legend('\it f\rm_z(z|d = +1)','\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('\it f\rm(z)');