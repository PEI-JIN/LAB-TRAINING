%%% BPSK Simulation Using MATLAB (Version 2)
%%% Assuming equal probability of the transmitted bits
%%% Written by P.-J. Su 2022/8/29

clc;    %% Clear command window
clear;  %% Remove items from workspace

EbOverN0_dB = -10:1:10;   %% Eb/N0 in dB
EbOverN0 = 10.^(EbOverN0_dB/10);    %% Eb/N0
EbOverN0_dB_approx = 0:1:8;    %% Eb/N0 in dB for approximation
EbOverN0_approx = 10.^(EbOverN0_dB_approx/10); %% Eb/N0 for approximation

Eb = 1;     %% Bit energy (J)
T = 1;      %% Bit interval (s)
num_bits = 10^7;    %% Number of data bit
N0 = Eb./EbOverN0;
sigma = sqrt(T*N0);
BPSK = [-sqrt(2*Eb*T),sqrt(2*Eb*T)];   %% Baseband modulator
P1 = 0.5;   %% The a priori probability: P(d = +1) = P(d = -1) = 0.5
r_th = 0.5*log((1 - P1)/P1)*(T*N0)/sqrt(2*Eb*T);   %% The threshold value
num_bit_error = zeros(size(EbOverN0_dB));     %% Number of error bits

% The theoretical BER of the optimal detector
BER_theory = qfunc(sqrt(2*EbOverN0));

% Approximation (8.16)
BER_approx = 0.5*exp(-EbOverN0_approx)./sqrt(pi*EbOverN0_approx);

for snrIdx = 1: length(EbOverN0_dB)
    tx_bit = rand(1, num_bits) > P1;  %% Transmitted data bit
    tx_sym = BPSK(tx_bit(1:num_bits)+1);  %% Transmitted symbol
    
    %  AWGN channel
    rx_sym = tx_sym + randn(1, length(tx_sym))*sigma(snrIdx);  
    
    rx_bit = zeros(size(tx_bit)); %% Reset rx bit decisions 
    rx_bit(rx_sym > r_th(snrIdx)) = 1;   %% Decision of bit '1'
    err_pat = xor(tx_bit, rx_bit); %% Compare tx and rx bits
    
    % Error counting
    num_bit_error(snrIdx) = num_bit_error(snrIdx) + sum(err_pat);   
end

% The simulated BER of the optimal detector
BER_sim = num_bit_error/num_bits;

% Plot the error rate curves
clf;
semilogy(EbOverN0_dB, BER_theory, 'b');
hold on;
semilogy(EbOverN0_dB_approx, BER_approx, 'r');
hold on;
semilogy(EbOverN0_dB, BER_sim, '--om');
axis([-10,10,1e-4,1]);

legend('Actual','Approximation (8.16)','Simulation');
xlabel('10 log_{10}z');
ylabel('\it P_E \rm \rightarrow');
grid;