% BER performance of an OFDM system with BPSK modulation in a static
% frequency-selective channel, written by SU PEI-JIN 2023/3/9

% Clear Command Window
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Close all figures whose handles are visible.
close all
% Start stopwatch timer.
tic

num_bit = 32*10^6;  % Number of data bits

Eb = 1;     % Energy per symbol
BPSK = [-sqrt(Eb), sqrt(Eb)];   % BPSK symbol mapping

SNR_in_dB = 0:1:20;     % SNR (Eb/N0) in dB
SNR = 10.^(SNR_in_dB/10);     % SNR (Eb/N0)
num_bit_error = zeros(size(SNR_in_dB));    % Number of bit error
BER_theo = zeros(size(SNR_in_dB));

% Initialize the random number generator
rng('default');

% Bit generation
tx_bit = ceil(2.*rand(1, num_bit))-1;

% BPSK Mapper
tx_BPSK_sym = BPSK(tx_bit(1:num_bit)+1); %% Transmit signal

% Serial to Parallel Converter
FFT_size = 64;  % Number of subcarriers
num_OFDM_sym = num_bit/FFT_size;
tx_sig = reshape(tx_BPSK_sym,[FFT_size,num_OFDM_sym]);

% IFFT
tx_IFFT_sig = ifft(tx_sig)*sqrt(FFT_size);

% CP Insertion
CP_size = FFT_size/4;
tx_OFDM_sym = [tx_IFFT_sig((FFT_size-CP_size+1):end,:); tx_IFFT_sig];

% Parallel to Serial Converter
tx_OFDM_sig = reshape(tx_OFDM_sym,1,(FFT_size+CP_size)*num_OFDM_sym);

% Frequency-selective Channel (LTI Channel)
num_taps = 8;   % Number of discrete path delay
h = exp(-0.25*(0:num_taps-1));

H = fft(h,FFT_size);    
fre_sel_ch_sig = conv(tx_OFDM_sig,h);

for snr_Idx = 1:length(SNR_in_dB)
    % AWGN Channel
    sigma = sqrt(Eb/(2*SNR(snr_Idx)));
    noise = randn(1, length(fre_sel_ch_sig))*sigma + 1i*randn(1, length(fre_sel_ch_sig))*sigma;
    AWGN_ch_sig = fre_sel_ch_sig + noise;
    
    rx_sig = AWGN_ch_sig(1:(FFT_size+CP_size)*num_OFDM_sym);
    
    % Serial to Parallel Converter
    rx_OFDM_sig = reshape(rx_sig,[FFT_size+CP_size,num_OFDM_sym]);
    
    % CP Removal
    shift_FFT_window = 0;   % different starting time of FFT window
    rx_OFDM_sym = rx_OFDM_sig((CP_size + 1)-shift_FFT_window:end-shift_FFT_window,:);
    rx_OFDM_sym = circshift(rx_OFDM_sym,-shift_FFT_window,1);
    
    % FFT
    rx_FFT_sig = fft(rx_OFDM_sym)/sqrt(FFT_size);

    % Frequency domain equalizaion (FEQ)
    H_transpose = H.';
    rx_FEQ_sig = rx_FFT_sig(:,1:end)./H_transpose;
    
    % Parallel to Serial Converter
    rx_BPSK_sym = reshape(rx_FEQ_sig,1,FFT_size*num_OFDM_sym);
    
    % BPSK Demapper
    rx_bit=zeros(1,FFT_size*num_OFDM_sym);
    rx_bit(rx_BPSK_sym>0) = 1; % decision of bit '1'

    % BER analysis
    err_pat = xor(tx_bit, rx_bit);      % compare tx and rx bits
    num_bit_error(snr_Idx) = num_bit_error(snr_Idx) + sum(err_pat); % error counting
    
    % Theoretical BER curve
    BER_theo(snr_Idx) = mean(qfunc(abs(H).*sqrt(2*SNR(snr_Idx))));
end

% BER performance via simulation
BER_sim = num_bit_error/num_bit;    % Simulated BER

semilogy(SNR_in_dB, BER_theo,'r-o');
hold on;
semilogy(SNR_in_dB, BER_sim,'b--x');
axis([0 15 1e-6 1]);
title('BER curve of an OFDM system in frequency-selective channel');
legend('Theory (BPSK modulation)','Simulation (BPSK modulation)');
xlabel('\itE_{b}\rm/\itN\rm_{0} (dB)');
ylabel('BER');
grid;

% Read elapsed time from stopwatch.
toc

% Plot the impulse response of the LTI channel
figure;
stem((0:num_taps-1),h);
axis([-1 num_taps 0 1.2]);
title('LTI channel (Frequency-selective channel)');
xlabel('n');
ylabel('h[n]');

% Plot the frequency response of the LTI channel
figure;
norm_freq = (0:(FFT_size-1))/FFT_size - 0.5;
plot(norm_freq,real(fftshift(H)));
title('Frequency response of the LTI channel');
xlabel('\itf_{N} \rm(normalize frequency)');
ylabel('Magnitude');