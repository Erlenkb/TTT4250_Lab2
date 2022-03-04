%% Pseudo code for lab 5 - Acoustic Measurement Terchnique
 
% Author: Robin Andre Rørstadbotnen, 27.04.2020
% Modified, 23.02.2022
% DO NOT DISTRIBUTE.
 
%% PART TWO FREE-FIELD

% All equation are taken from lecutre notes, some input equations taken
% elsewhere. Note that Figure 6 in the method document should be improved
% as the heights should be more intuitive.

%% Input
T = 20;                              % Celcius
rho =   ;                            % kg/m3, density of air
c = 343*sqrt( ( T + 273 ) / 293 );   % m/s, sound speed in air


dd  = ;                              % Horizontal distance mic - loudpeaker
hh  = ;                              % Height loudspeaker
hm1 = ;                              % Height mic 1
hm2 = ;                              % Height mic 2

% From Figure 3 (lecture notes) ir Figure 2 (exercise sheet)
rd1 = sqrt ( (hh - hm1)^2 + dd^2 );  % Total direct distance travelled to mic 1
rd2 = sqrt ( (hh - hm2)^2 + dd^2 ); 
rr1 =                                % Total reflected distance travelled to mic 1
rr2 = 
theta =                              % Reflection angle

%% Importing measurements
% Quasi-Anechoic measurements

tmax      =                        % Extract only direct and the primary reflected puls from the rockwool
path2data =                        % Path to the folder containing the data
FFfile1   = '.etx';                % Or .txt  
FFfile2   = '.etx';                % Or .txt

Psamp1_in = importdata([path2data FFfile1],'\t',22);
Psamp2_in = importdata([path2data FFfile2],'\t',22); 

tt        = Psamp1_in(:,1) ;
dt        = tt(2) - tt(1)  ;       % Sample time
fs        = 1/dt;                  % Sample frequency

idx_tmax  = tmax/dt;               % or, find( tt == tmax ) ;

% Extract signal of interest
p1 = Psamp1_in.data(1:idx_tmax,2);   
p2 = ;


% Quality control of input signal, also for prep questions
figure
plot()
ylabel()
xlabel()
legend()
grid on
title()
set(gca,'FontSize',12,'Fontweight','bold')

%% Frequency analysis, also for prep
n  = ;
ff = ;
P1 = fft(p1, ) / ( ) ;
P2 = fft(p2, ) / ( ) ;

% Plot it in dB
figure 
semilogx()
xlabel 
set(gca,'FontSize',20,'Fontweight','bold') % ++++

%% Computation

% Transfer function - eq.25 (remember component wise mulp)

H12 = ;

k = ;           % Wave number

R = ;

Z = ;

alpha = ;

% Rest (reflection coefficient, acoustic impedance) from 6.2, absorption
% coefficient from 3
% Be carefult that you have all the parameters for the reflection
% coefficient in the right place and enough parentheses

%% Plot the result


%% Save results


















    
