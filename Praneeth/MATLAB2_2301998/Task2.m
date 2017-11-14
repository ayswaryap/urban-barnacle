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

%% Question 1  

% Assume WGN noise in the measurements
%%%%% Here is the hypothesis model:- H0: a = n and H1: a = X + n
%%% We detect the noise using rayleigh distribution 


% Unmwanted data filtering  
filterData = SigPwr(find(SigPwr~=-Inf)); 

%Average Noise Power 
noiseFloor = mean(filterData); 

%Noise Variance
noiseVar = noiseFloor - min(filterData);      

%Threshold - calculated as the sum of noise Floor and the noise standard
%deviation
Threshold = noiseFloor + noiseVar;       

% Probability of False Alarm(PFA) - when noise voltage exceeds the threshold voltage 
PFA = exp(-0.5*10^(0.1 * (Threshold - noiseFloor) ) );

% To Print the Probability of False Alarm, Threshold and Noise Floor
fprintf('PFA=  %f ,Threshold =  %f dBm, Noise Floor  = %f.\n',PFA,Threshold,noiseFloor)


%%
%%Question 2

H0 = find(SigPwr<Threshold);
H1 = find(SigPwr>=Threshold);
Sample(H0)=0;
Sample(H1)=1;
figure(2);

subplot(3,1,1) % To break the figure into 3 x 1 matrix of small axes and 1 gives the current plot H0
plot(SigPwr(H0));
xlabel('Time(sec)');
ylabel('Signal Power (dBm)');
title('H0');

hold all;

subplot(3,1,2) % To break the figure into 3 x 1 matrix of small axes and 2 gives the current plot H1
plot(SigPwr(H1));
xlabel('Time(sec)');
ylabel('Signal Power (dBm)');
title('H1');

subplot(3,1,3) % To break the figure into 3 x 1 matrix of small axes and 3 gives the current plot Sample
plot(Sample); %Gives the square signal
xlabel('Time(sec)');
ylabel('Signal Power (dBm)');
title('Square signal');



%%
%Question 3 and 4


% clear all;
% clc;
% 
% % Loading the given data file
% load('data_high_snr.mat');
% 
% figure(1);
% %Store the given data_high_snr
% SigPwr = 10*log10(abs(data).^2); 
% 
% %Plot
% plot(timet,SigPwr);
% xlabel('Time(sec)');
% ylabel('Signal Power (dBm)');
% hold all;

% Differentiator (filter) to find pulse duration and periodicity
diff_signal = diff(data);
diff_signal_sq = abs(diff_signal).^2;

mean_diff_sig_sq = mean(diff_signal_sq);
std_diff_sig_sq = std(diff_signal_sq);

% threshold constant (32)
th_factor = 32;
th_value = mean_diff_sig_sq + std_diff_sig_sq * th_factor; % threshold helps us in reducing the noise variance

th_value_dBm = 10*log10(th_value);

% finding peak locations based on threshold (includes both raising and
% falling edge (odd is raising edge, even is falling edge)
peak_locations = find(diff_signal_sq >= th_value);

% Using diff will skip one sample, therefore compensating it by adding 1
raise_locs = timet(peak_locations(1:2:end - 1) + 1);
fall_locs = timet(peak_locations(2:2:end) - 1);

signal_duration = mean(fall_locs - raise_locs);
signal_periodicity = mean([diff(fall_locs), diff(raise_locs)]);

fprintf('Pulse Duration=  %f sec,Pulse Period =  %f sec, threshold = %fdB.\n',signal_duration,signal_periodicity, th_value_dBm)

%% Question 5 
%Detecting low SINR (weak signals)

figure(3);
% 0.001 frame boundary (accumulating over 16 frames)
rsh_data = reshape(data,16,[]);

subplot(2,1,1);plot(timet(1:length(rsh_data(1,:))),db(abs(rsh_data(1,:)),'voltage'));
subplot(2,1,1);xlabel('Time Duration in Sec');ylabel('Signal power in dBm');title('Low SNR behavior');

% non coherent addition (to avoid deterioration due to frequency offset
% across multple frames)

avg_signal = mean(abs(rsh_data),1);
subplot(2,1,2);plot(timet(1:length(rsh_data(1,:))),db(avg_signal,'voltage'));
subplot(2,1,2);xlabel('Time Duration in Sec');ylabel('Signal power in dBm');title('Noise supression by noncoherent combining');





















