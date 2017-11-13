clear all
clc
load('data_high_snr.mat')     %%%%%% Loading data file
plot(timet,10*log10(abs(data).^2))
%%% Assume WGN noise in the measurements
%%%%% Here is the hypothesis model
% H0: r = w
% H1: r = A + w
%%% Then using Neyman-Pearson Theorem %%%%
%%%  P(r|H1)/P(r|H0) > eta
%%% for simplicity we use -> eta = Pr(H0)/Pr(H1) 
%%% Then we detect H1 if -> r > N/A*log(eta) + A/2 
%%% We can convert this to dBm ---> 10*log10(abs(r).^2)
%%% No of pulse samples around 1300 -> Pr(H1) = 0.0065 
%%% From the first plot --> N = 10^9 mW ; A = 10^3 mV ? 
%%% Therefore we use our threshold -66dBm

