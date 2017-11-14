% ??????????????????????????????????????????????????????????????????????
%          STATISTICAL SIGNAL PROCESSING (Matlab Exercise #2)
%                     (Least Square estimation)
%
% ??????????????????????????????????????????????????????????????????????

% ***NOTE***: use the button "Run Section" to visualize the diagram of the
% recorded positions of 'a' with respect to 'B'.



% (X,Y) Measurements of the position of the celestial body 'a' for each
% month:
X = [69.9610   69.4111   68.1078   68.4906   66.5073   65.0106   63.3201   60.0063   57.8493   54.9152   53.1017   49.5782];
Y = [20.7559   22.6554   24.3398   24.9237   27.6180   29.2726   29.8596   31.5923   32.9917   35.2425   36.1634   36.1643];

% Orbital period of 'a':
T = 12*6;

% Plot the diagram:
scatter(X, Y, '.', 'LineWidth',1);
hold on;
h = ellipse(5, 5, 0, 0,0, 'r');
set(h, 'LineWidth', 2,  'LineStyle','-');
hold off;
axis equal;
xlabel X;
ylabel Y;
text(5,-5, 'B', 'Color', 'Red', 'FontWeight', 'Bold', 'FontSize', 16);
text(X(1)-2,Y(1)-1, 'a', 'Color', 'Blue', 'FontWeight', 'Bold', 'FontSize', 16);


%% 
% ***********************************************************************************
% Piece of code that estimates the four parameters of the ellipse using LS-estimation
% ***********************************************************************************

% ** LS Linear Estimator :- Y = HX+N is X^{hat} = inv(H'*H)*H'*Y ****

%************ Observation matrices******************************

k = 12; %Given
Obs1 =[]; %H0
Obs2 =[]; %H1

for K = 1:k % recodered in 12 months 
    Obs1 = [Obs1; cos(2*pi*K/T) 1]; %T is the orbital period of a defined earlier
    Obs2 = [Obs2; sin(2*pi*K/T) 1];
end



%******** LS estimate: R1,C1 and R2,C2**********************************

LS_R1 = inv(Obs1' * Obs1) * Obs1' * X'; %LS estimate for X_k
LS_R2 = inv(Obs2' * Obs2) * Obs2' * Y'; %LS estimate of Y_k


%  NOTE: The estimated parameters R1,R2,C1,C2 must be stored in a 4-element vector variable named P.


P = [LS_R1(1) LS_R2(1) LS_R1(2) LS_R2(2)]

h = ellipse(P(1), P(2), 0, P(3), P(4), 'b');
set(h, 'LineStyle','--');

