
clear all;
clc;

% Loading the given data file
load('data_high_snr.mat');

figure(1);
%Store the given data_high_snr
SigPwr = 10*log10(abs(data).^2); 

%Plot
plot(timet,SigPwr);
xlabel('Time(sec)');
ylabel('Signal Power (dBm)');
hold all;

% Differentiator (filter) to find pulse duration and periodicity
diff_signal = diff(data);
diff_signal_sq = abs(diff_signal).^2;

mean_diff_sig_sq = mean(diff_signal_sq);
std_diff_sig_sq = std(diff_signal_sq);

% threshold constant (32)
th_factor = 32;
th_value = mean_diff_sig_sq + std_diff_sig_sq * th_factor;

% finding peak locations based on threshold (includes both raising and
% falling edge (odd is raising edge, even is falling edge)
peak_locations = find(diff_signal_sq >= th_value);

% Using diff will skip one sample, therefore compensating it by adding 1
raise_locs = timet(peak_locations(1:2:end - 1) + 1);
fall_locs = timet(peak_locations(2:2:end) - 1);

signal_duration = mean(fall_locs - raise_locs);
signal_periodicity = mean([diff(fall_locs), diff(raise_locs)]);

fprintf('Pulse Duration=  %f sec,Pulse Period =  %f sec.\n',signal_duration,signal_periodicity)

%% Detecting low SINR (weak signals)

figure(2);
% 0.001 frame boundary (accumulating over 16 frames)
rsh_data = reshape(data,16,[]);

subplot(2,1,1);plot(timet(1:length(rsh_data(1,:))),db(abs(rsh_data(1,:)),'voltage'));
subplot(2,1,1);xlabel('Time Duration in Sec');ylabel('Signal power in dBm');title('Low SNR behavior');

% non coherent addition (to avoid deterioration due to frequency offset
% across multple frames)

avg_signal = mean(abs(rsh_data),1);
subplot(2,1,2);plot(timet(1:length(rsh_data(1,:))),db(avg_signal,'voltage'));
subplot(2,1,2);xlabel('Time Duration in Sec');ylabel('Signal power in dBm');title('Noise supression by noncoherent combining');






