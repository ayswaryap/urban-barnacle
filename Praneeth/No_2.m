clear all
clc
load('data_high_snr.mat')                   %%%%%% Loading data file
figure(1)
plot(timet,10*log10(abs(data).^2))
X = 10*log10(abs(data).^2);                 %%%%%% Store data in dBm
Y = X(find(X~=-Inf));                       %%%%%% Filtering unwanted data.  
%% Assume WGN noise in the measurements
%%%%% Here is the hypothesis model
% H0: r = w
% H1: r = A + w
%%% Since noise is complex gausian we detect noise with rayleigh distribution 
%%% We have the derivation in the report for PFA
%%% Set some values;
N_Floor = mean(Y);              %%%%%%%%%%%% Noise floor ->> average noise power
N_var = N_Floor - min(Y);       %%%%%%%%%%%% Noise variance
P_Th = N_Floor + N_var;         %%%%%%%%%%%% Threshold

%%%%%% From theoretical formula - PFA
PFA = exp(-0.5*10^(0.1 * (P_Th - N_Floor) ) );
fprintf('PFA=  %f ,Threshold =  %f, Noise Floor  = %f.\n',PFA,P_Th,N_Floor);

%%%%%%%%%%% Second Part %%%%%%%%%%%%%%%%
H0 = find(X<P_Th);
H1 = find(X>=P_Th);
Z(H0)=0;
Z(H1)=1;
figure(2)
subplot(3,1,1)
plot(X(H0));
subplot(3,1,2)
plot(X(H1));
subplot(3,1,3)
plot(Z);
%%%%%%%%%%% Next Part %%%%%%%%%%%%%%%%%%%%
Fs = length(timet)/(timet(end)-timet(1));     % Sampling frequency
[P_val P_loc P_duration] = findpeaks(Z,Fs);   % find pulse duration and pulse frequency
%%%%% Remove false Peaks %%%%%%%%%
N_limit = 5/Fs;
j=0;
for i=1:length(P_duration)
    if P_duration(i) > N_limit
        j=j+1;
        Duration_vector(j)= P_duration(i);
        loc_vector(j) = P_loc(j);
    end
end
Pulse_Duration = mean(Duration_vector)
A=[];
for i=2:j
    A=[A loc_vector(j)-loc_vector(j-1)];
end
Pulse_Period = mean(A)
