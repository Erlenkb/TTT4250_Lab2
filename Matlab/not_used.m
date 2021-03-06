%% Pseudo code for lab 5 - Acoustic Measurement Terchnique
 
% Author: Robin Andre Rørstadbotnen, 27.04.2020
% Modified, 23.02.2022
% DO NOT DISTRIBUTE.
 
%% PART TWO FREE-FIELD

% All equation are taken from lecutre notes, some input equations taken
% elsewhere. Note that Figure 6 in the method document should be improved
% as the heights should be more intuitive.

%% Input
close all
clear 
clc
rho_0 = 1.186;
p_0 = 101.325;

T = 20;                              % Celcius
rho = 1.225;                            % kg/m3, density of air
c = 343.2*sqrt( ( T + 271 ) / 293 );   % m/s, sound speed in air

%1.004041 / c


dd  = 1;                              % Horizontal distance mic - loudpeaker
hh  = 0.39;                              % Height loudspeaker
hm1 = 0.45;                              % Height mic 1
hm2 = 0.5;                              % Height mic 2

% From Figure 3 (lecture notes) ir Figure 2 (exercise sheet)
rd1 = sqrt ( (hh - hm1)^2 + dd^2 );  % Total direct distance travelled to mic 1
rd2 = sqrt ( (hh - hm2)^2 + dd^2 ); 
rr1 = sqrt ( (hh + hm1)^2 + dd^2);   % Total reflected distance travelled to mic 1
rr2 = sqrt ( (hh + hm2)^2 + dd^2); 
rr  = sqrt ( (hh + ((hm1+hm2)/2))^2 + dd^2);
theta1 = (pi/2 - acos(1/rr1));                              % Reflection angle
theta2 = (pi/2 - acos(1/rr2));
theta  = (pi/2 - acos(1/rr));

%% Importing measuremen
% Quasi-Anechoic measurements

tmax      = 0.20;                        % Extract only direct and the primary reflected puls from the rockwool
%path2data =                        % Path to the folder containing the data
FFfile1   = 'Free_Field_45cm_Height_d1m.etx';                % Or .txt  
FFfile2   = 'Free_Field_50cm_Height_d1m.etx';                % Or .txt

Psamp1_in = importdata([FFfile1],'\t',22);
Psamp2_in = importdata([FFfile2],'\t',22); 

tt        = Psamp1_in.data(:,1) ;
dt        = tt(2) - tt(1);       % Sample time
fs        = 1/dt;                  % Sample frequency

idx_tmax  = 8e-3/(tt(2)-tt(1));               % or, find( tt == tmax ) ;

% Extract signal of interest
p1 = Psamp1_in.data(1:idx_tmax,2);   
p2 = Psamp2_in.data(1:idx_tmax,2);
tt = tt(1:idx_tmax);

% Quality control of input signal, also for prep questions
figure(32)
subplot(2,2,1)
plot(tt, p1)
hold on
plot(tt, p2)
hold off
ylabel("Magnitude")
xlabel("Time [ms]")
legend()
grid on
title("p1")
set(gca,'FontSize',12,'Fontweight','bold')
set(gcf,'units','centimeters','position',[2,1,29.7,21.0])




%% Frequency analysis, also for prep
n  = 2^nextpow2( size(p1,1) );  
paddingnumber = n - size(p1,1);
p1 = padarray(p1, paddingnumber, 0, "post");
p2 = padarray(p2, paddingnumber, 0, "post");

ff = fs*(0:(n-1))/n;
frecvec1 = fft(p1,n);
frecvec2 = fft(p2,n);

% Plot it in dB
figure 
semilogx(ff, 20*log10(frecvec1))
hold on
grid on
semilogx(ff, 20*log10(frecvec2))
hold off
xlabel("Frequency [Hz]");
xlim([100,2000])
set(gca,'fontsize',12,'fontweight','bold'); % ++++

%% Computation

% Transfer function - eq.25 (remember component wise mulp)

H12_Free = transpose(sqrt(( p2./p1).^2));
%figure(12)
%plot(ff, H12_Free)
%title("Transfer function")
%xlim([100 2000])


ww = 2*pi*ff;



k =  ww / c;           % Wave number

R_num = ((exp(-1i*k*rd2))/rd2)   -  H12_Free .*((exp(-1i*k*rd1))/rd1);
R_den = H12_Free.*((exp(-1i*k*rr1))/rr1)   - ((exp(-1i*k*rr2))/rr2);

R = R_num ./ R_den;
figure(32)
subplot(2,2,2)
semilogx(ff, abs(R))
xlim([100 2000])
ylim([0 1])
title("Reflection")
grid on


%% Stegvis utregning av Impedansen med og uten skalaren
scalar = (-1*rho*c) / cos(theta);
Z = (R+1) ./ (R-1);
Z = scalar * Z;

%% Formel for Impedanse
%Z = (-1*c*rho*(1+R)) ./ (cos(theta)*(R-1));


%Z = -1*(rho*c*(R+1)*sec(theta)) ./ (R-1);
alpha = 1 - abs(R).^2;
figure(32)
subplot(2,2,3)
semilogx(ff, alpha)
xlim([100 2000])
ylim([0 1])
grid on
title("Absoprtion")

figure(32)
subplot(2,2,4)
semilogx(ff,real(Z),'r')
hold on
semilogx(ff,imag(Z),'b');
semilogx(ff,abs(Z),'k');
grid on
% ylim([-4 4])
xlim([100 2000]);
title('Impedance Z')
xlabel('Frequency [Hz]' )
legend('')
hold off
exportgraphics(figure(32), ['Free_Field.png'],'Resolution',450)

% Rest (reflection coefficient, acoustic impedance) from 6.2, absorption
% coefficient from 3
% Be carefult that you have all the parameters for the reflection
% coefficient in the right place and enough parentheses

%% Plot the result


%% Save results


















    
