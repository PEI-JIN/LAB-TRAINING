%%% Scatter Plot and Histogram @ SNR = 3dB and SNR = 10dB
%%% Written by P.-J. Su 2022/9/5

clc;    %% Clear command window
clear;  %% Remove items from workspace

Eb = 1;     %% Bit energy (J)
T = 1;      %% Bit interval (s)
num_bits = 10^6;    %% Number of data bit
P1 = 0.5;   %% The a priori probability: P(1) = 0.5

%%% Transmitted data bit %%%
tx_bit = rand(1, num_bits) > (1-P1);

%%% Transmitted symbol %%%
bit_1_Idx = 1;
bit_0_Idx = 1;
tx_1_sym = zeros(1,sum(tx_bit) + 1);
tx_0_sym = zeros(1,num_bits - sum(tx_bit) + 1);
for bitIdx = 1:num_bits
    if tx_bit(1, bitIdx)
        tx_1_sym(bit_1_Idx) = sqrt(2*Eb*T);
        bit_1_Idx = bit_1_Idx+1;
    else
        tx_0_sym(bit_0_Idx) = -sqrt(2*Eb*T);
        bit_0_Idx = bit_0_Idx+1;
    end        
end
    
%%%  AWGN channel (SNR = 3dB) %%%
SNR_3dB = 10.^(3/10);
N0_3dB = 2*Eb/SNR_3dB;  %% SNR = 2*Eb/N0
sigma_3dB = sqrt(N0_3dB*T);

% When data bit d = +1 is transmitted (Source bit: 1)
z_1_3dB = tx_1_sym + randn(size(tx_1_sym))*sigma_3dB;
% When data bit d = -1 is transmitted (Source bit: 0)
z_0_3dB = tx_0_sym + randn(size(tx_0_sym))*sigma_3dB;

%%%  AWGN channel (SNR = 10dB) %%%
SNR_10dB = 10.^(10/10);
N0_10dB = 2*Eb/SNR_10dB;    %% SNR = 2*Eb/N0
sigma_10dB = sqrt(N0_10dB*T);

% When data bit d = +1 is transmitted (Source bit: 1)
z_1_10dB = tx_1_sym + randn(size(tx_1_sym))*sigma_10dB;
% When data bit d = -1 is transmitted (Source bit: 0)
z_0_10dB = tx_0_sym + randn(size(tx_0_sym))*sigma_10dB;

%%% Plot the scatter plot and histogram @ SNR = 3dB %%%
clf;
tiledlayout(2,1);
% Top plot
ax1 = nexttile;
scatter(ax1,z_1_3dB,zeros(size(tx_1_sym)),'blue','filled','o');
xlim([-7 7]);
title('Scatter Plot @ SNR = 3dB (d = +1)');
yticks([]);
xlabel('z');
% Bottom plot
ax2 = nexttile;
scatter(ax2,z_0_3dB,zeros(size(tx_0_sym)),'red','filled','d');
xlim([-7 7]);
title('Scatter Plot @ SNR = 3dB (d = -1)');
yticks([]);
xlabel('z');

figure;
nbins = 30;
h1 = histogram(z_1_3dB,nbins);
hold on;
h2 = histogram(z_0_3dB,nbins);
title('Histogram @ SNR = 3dB (Count)');
legend('d = +1','d = -1');
xlabel('z');
ylabel('Count');

d_1_data_3dB = ([h1.BinEdges 0]+[0 h1.BinEdges])/2;
d_1_x_3dB = d_1_data_3dB(2:nbins+1);
d_1_count_3dB = h1.Values;
d_1_dx_3dB = h1.BinWidth;
d_1_prob_3dB = d_1_count_3dB/num_bits;
d_1_pdf_3dB = d_1_prob_3dB/d_1_dx_3dB;

d_0_data_3dB = ([h2.BinEdges 0]+[0 h2.BinEdges])/2;
d_0_x_3dB = d_0_data_3dB(2:nbins+1);
d_0_count_3dB = h2.Values;
d_0_dx_3dB = h2.BinWidth;
d_0_prob_3dB = d_0_count_3dB/num_bits;
d_0_pdf_3dB = d_0_prob_3dB/d_0_dx_3dB;

%%% Theory %%%
x_1 = (sqrt(2*Eb*T)-3):.1:(sqrt(2*Eb*T)+3);
x_0 = (-sqrt(2*Eb*T)-3):.1:(-sqrt(2*Eb*T)+3);

d_1_theo_3dB = P1*(1/(sqrt(2*pi)*sigma_3dB))* ...
               exp((-0.5/N0_3dB)*(x_1-sqrt(2*Eb*T)).*(x_1-sqrt(2*Eb*T)));
d_0_theo_3dB = (1-P1)*(1/(sqrt(2*pi)*sigma_3dB))* ...
               exp((-0.5/N0_3dB)*(x_0+sqrt(2*Eb*T)).*(x_0+sqrt(2*Eb*T)));

figure;
bar(d_1_x_3dB,d_1_prob_3dB);
hold on;
bar(d_0_x_3dB,d_0_prob_3dB);
ylim([0 0.1]);
title('Histogram @ SNR = 3dB (Probability)');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('Probability');

figure;
bar(d_1_x_3dB,d_1_pdf_3dB);
hold on;
bar(d_0_x_3dB,d_0_pdf_3dB);
hold on;
plot(x_1, d_1_theo_3dB,'--b');
hold on;
plot(x_0, d_0_theo_3dB,'--r');
ylim([0 1]);
title('Histogram @ SNR = 3dB (pdf)');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1) \bf (Simulation)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1) \bf (Simulation)', ...
       '\it P\rm_1\cdot\it f\rm_z(z|d = +1) \bf (Theory)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1) \bf (Theory)');
xlabel('z');
ylabel('\it f\rm(z)');

figure;
plot(x_1, d_1_theo_3dB,'--b');
hold on;
plot(x_0, d_0_theo_3dB,'--r');
ylim([0 1]);
title('Theoretical behavior of the a posteriori probability @ SNR = 3dB');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('\it f\rm(z)');

%%% Plot the scatter plot and histogram @ SNR = 10dB %%%
figure;
tiledlayout(2,1);
% Top plot
ax1 = nexttile;
scatter(ax1,z_1_10dB,zeros(size(tx_1_sym)),'blue','filled','o');
xlim([-4 4]);
title('Scatter Plot @ SNR = 10dB (d = +1)');
yticks([]);
xlabel('z');
% Bottom plot
ax2 = nexttile;
scatter(ax2,z_0_10dB,zeros(size(tx_0_sym)),'red','filled','d');
xlim([-4 4]);
title('Scatter Plot @ SNR = 10dB (d = -1)');
yticks([]);
xlabel('z');

figure;
h3 = histogram(z_1_10dB,nbins);
hold on;
h4 = histogram(z_0_10dB,nbins);
title('Histogram @ SNR = 10dB (Count)');
legend('d = +1','d = -1');
xlabel('z');
ylabel('Count');

d_1_data_10dB = ([h3.BinEdges 0]+[0 h3.BinEdges])/2;
d_1_x_10dB = d_1_data_10dB(2:nbins+1);
d_1_count_10dB = h3.Values;
d_1_dx_10dB = h3.BinWidth;
d_1_prob_10dB = d_1_count_10dB/num_bits;
d_1_pdf_10dB = d_1_prob_10dB/d_1_dx_10dB;

d_0_data_10dB = ([h4.BinEdges 0]+[0 h4.BinEdges])/2;
d_0_x_10dB = d_0_data_10dB(2:nbins+1);
d_0_count_10dB = h4.Values;
d_0_dx_10dB = h4.BinWidth;
d_0_prob_10dB = d_0_count_10dB/num_bits;
d_0_pdf_10dB = d_0_prob_10dB/d_0_dx_10dB;

%%% Theory %%%
d_1_theo_10dB = P1*(1/(sqrt(2*pi)*sigma_10dB))* ...
               exp((-0.5/N0_10dB)*(x_1-sqrt(2*Eb*T)).*(x_1-sqrt(2*Eb*T)));
d_0_theo_10dB = (1-P1)*(1/(sqrt(2*pi)*sigma_10dB))* ...
               exp((-0.5/N0_10dB)*(x_0+sqrt(2*Eb*T)).*(x_0+sqrt(2*Eb*T)));

figure;
bar(d_1_x_10dB,d_1_prob_10dB);
hold on;
bar(d_0_x_10dB,d_0_prob_10dB);
ylim([0 0.1]);
title('Histogram @ SNR = 10dB (Probability)');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('Probability');

figure;
bar(d_1_x_10dB,d_1_pdf_10dB);
hold on;
bar(d_0_x_10dB,d_0_pdf_10dB);
hold on;
plot(x_1, d_1_theo_10dB,'--b');
hold on;
plot(x_0, d_0_theo_10dB,'--r');
ylim([0 1]);
title('Histogram @ SNR = 10dB (pdf)');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1) \bf (Simulation)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1) \bf (Simulation)', ...
       '\it P\rm_1\cdot\it f\rm_z(z|d = +1) \bf (Theory)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1) \bf (Theory)');
xlabel('z');
ylabel('\it f\rm(z)');

figure;
plot(x_1, d_1_theo_10dB,'--b');
hold on;
plot(x_0, d_0_theo_10dB,'--r');
ylim([0 1]);
title('Theoretical behavior of the a posteriori probability @ SNR = 10dB');
legend('\it P\rm_1\cdot\it f\rm_z(z|d = +1)', ...
       '(1-\it P\rm_1)\cdot\it f\rm_z(z|d = -1)');
xlabel('z');
ylabel('\it f\rm(z)');

