%%% Appendix 1 Simulation Program
%%% Written by P.-J. Su 2022/8/24

% Clear Command Window.
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Start stopwatch timer.
tic

SNR_dB = 0:2:14;   %% SNR in dB
SNR = 10.^(SNR_dB/10);

Pe_theoretical  = qfunc(sqrt(SNR));
Pe_approximation_1 = 0.5*exp(-0.5*sqrt(SNR).*sqrt(SNR));
constant_1 = 1/sqrt(2*pi);
Pe_approximation_2 = constant_1*exp(-0.5*sqrt(SNR).*sqrt(SNR))./sqrt(SNR);

clf;
semilogy(SNR_dB, Pe_theoretical, 'b', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, Pe_approximation_1, 'Color', '#008000', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, Pe_approximation_2, 'r', 'LineWidth', 2);
axis([0,14,1e-7,1]);

legend('Theory (ML Detector)','Approximation 1','Approximation 2');
xlabel('SNR (dB)');
ylabel('BER');
grid;

% Read elapsed time from stopwatch.
toc