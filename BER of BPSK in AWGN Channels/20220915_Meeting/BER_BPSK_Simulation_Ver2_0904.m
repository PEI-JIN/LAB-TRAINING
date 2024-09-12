%%% BPSK Simulation Using MATLAB (Version 2)
%%% Written by P.-J. Su 2022/9/4

clc;    %% Clear command window
clear;  %% Remove items from workspace

SNR_dB = 0:2:14;    %% SNR (= 2Eb/N0) in dB
EbOverN0_dB = SNR_dB - 10*log10(2);   %% Eb/N0 in dB
EbOverN0 = 10.^(EbOverN0_dB/10);    %% Eb/N0
EbOverN0_dB_approx = 0:1:8;    %% Eb/N0 in dB for approximation
EbOverN0_approx = 10.^(EbOverN0_dB_approx/10); %% Eb/N0 for approximation

Eb = 1;     %% Bit energy (J)
T = 1;      %% Bit interval (s)
num_bits = 10^8;    %% Number of data bits
N0 = Eb./EbOverN0;
sigma = sqrt(T*N0);
BPSK = [-sqrt(2*Eb*T),sqrt(2*Eb*T)];   %% Baseband modulator
P1 = 0.5;   %% The a priori probability: P(1) = 0.5
r_th = 0.5*log((1 - P1)/P1)*(T*N0)/sqrt(2*Eb*T);   %% The threshold value
num_error_MAP = zeros(size(EbOverN0_dB));   %% Number of error bits (MAP)
num_error_ML = zeros(size(EbOverN0_dB));   %% Number of error bits (ML)
num_tx_bits = 10^5;     %% Number of tx bits

% The theoretical BER of the MAP detector
MAP_constant_1 = sqrt(2*Eb*T) - r_th;
MAP_constant_2 = r_th + sqrt(2*Eb*T);
BER_theory_MAP = P1*qfunc(MAP_constant_1./sigma) + ...
                 (1-P1)*qfunc(MAP_constant_2./sigma);

% The theoretical BER of the ML detector
BER_theory_ML = qfunc(sqrt(2*EbOverN0));

% Approximation (8.16)
BER_approx = 0.5*exp(-EbOverN0_approx)./sqrt(pi*EbOverN0_approx);

for numTest = 1: num_bits/num_tx_bits
    fprintf(1,'Number of test = %d \n', numTest);
    for snrIdx = 1: length(EbOverN0_dB)        
        tx_bit = rand(1, num_tx_bits) > (1-P1);  %% Transmitted data bit
        tx_sym = BPSK(tx_bit(1:num_tx_bits) + 1);  %% Transmitted symbol
        
        %  AWGN channel
        rx_sym = tx_sym + randn(1, length(tx_sym))*sigma(snrIdx);  
        
        %%% MAP Detector %%%
        rx_bit_MAP = zeros(size(tx_bit));   %% Reset rx MAP decisions 
        rx_bit_MAP(rx_sym > r_th(snrIdx)) = 1;   %% MAP Decision of bit '1'
        % Compare tx and rx bits (MAP)
        err_pat_MAP = xor(tx_bit, rx_bit_MAP); 
    
        % Error counting (MAP criterion)
        num_error_MAP(snrIdx) = num_error_MAP(snrIdx) + sum(err_pat_MAP);
    
        %%% ML Detector %%%
        rx_bit_ML = zeros(size(tx_bit));   %% Reset rx ML decisions 
        rx_bit_ML(rx_sym > 0) = 1;   %% ML Decision of bit '1'
        % Compare tx and rx bits (ML)
        err_pat_ML = xor(tx_bit, rx_bit_ML);  
    
        % Error counting (ML criterion)
        num_error_ML(snrIdx) = num_error_ML(snrIdx) + sum(err_pat_ML);
    end
end
fprintf(1,'Simulation Completed\n');

% The simulated BER of the optimal detector (i.e., MAP)
BER_sim_MAP = num_error_MAP/num_bits;

% The simulated BER of the ML detector
BER_sim_ML = num_error_ML/num_bits;

% Plot the error rate curves
clf;
semilogy(EbOverN0_dB, BER_theory_MAP, '-^m');
hold on;
semilogy(EbOverN0_dB, BER_sim_MAP, '--ob');
hold on;
semilogy(EbOverN0_dB, BER_theory_ML, '-*g');
hold on;
semilogy(EbOverN0_dB, BER_sim_ML, '--sr');
hold on;
semilogy(EbOverN0_dB_approx, BER_approx, 'k');
x_axis_min = 0 - 10*log10(2);
x_axis_max = 14 - 10*log10(2);
axis([x_axis_min,x_axis_max,1e-7,5]);

legend('Theory (MAP Detector)','Simulation (MAP Detector)', ...
       'Theory (ML Detector)','Simulation (ML Detector)', ...
       'Approximation (8.16)');
xlabel('E_b/N_0 (dB)');
ylabel('BER');
grid;