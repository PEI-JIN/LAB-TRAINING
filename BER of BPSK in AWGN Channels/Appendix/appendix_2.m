%%% Appendix 2 Simulation Program
%%% Written by P.-J. Su 2022/8/24

% Clear Command Window.
clc
% Remove all variables from the current workspace,
% releasing them from system memory.
clear
% Start stopwatch timer.
tic

EbOverN0_dB = -10:1:10;   %% Eb/N0 in dB
EbOverN0 = 10.^(EbOverN0_dB/10);

Pe_theoretical  = qfunc(sqrt(2*EbOverN0));
Pe_approximation(1,11:21) = 0.5*exp(-EbOverN0(1,11:21))./sqrt(pi*EbOverN0(1,11:21));

clf;
semilogy(EbOverN0_dB, Pe_theoretical, EbOverN0_dB, Pe_approximation);

axis([-10,10,1e-4,1]);
text(-6, 0.3,'Actual');
text(3, 0.04,'Approximation (8.16)');

xlabel('10 log_{10}z');
ylabel('P_E \rightarrow');

% Read elapsed time from stopwatch.
toc