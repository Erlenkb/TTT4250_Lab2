%% Pseudo code for lab 5 - Acoustic Measurement Terchnique
 
% Author: Robin Andre Rørstadbotnen, 27.04.2020
% Modified, 23.02.2022
% DO NOT DISTRIBUTE.

close all
clear 
clc
 
% PART ONE IMPEDANCE TUBE 
% All equations from ISO10534-2
% Input (section 7.2 and method document):
T   = 291;     % Assumed temperature
rho = 1.225;
c   = 343.2*sqrt(T/293);       % eq. (5)
 
x_1 = 0.180;       %Distance mic 1
d   = 0.1;         %Diameter of tube
s   = 0.08;        %Distance between microphones
 
% Working frequency range 4.2
% f_l < f < f_u
f_l =  100;       % .05*c/s;   
f_u = .58*c/d;     % Eq. 2 should this be eq.3?
if .45*c/s <= f_u
    f_u = .58*c/d;         % Eq. 3
end
 
% Load .etx files to Matlab
path2files = "C:\Users\erlen\Downloads\Lab2_TTT4250"
fileName1  = "Imp_tube_12.etx"
fileName2  = "Imp_tube_21.etx"

Tube_f1    = fileName1;
Tube_f2    = fileName2;
Psamp1     = importdata(Tube_f1,'\t',22); % 22 is the number of header lines
Psamp2     = importdata(Tube_f2,'\t',22);

% Extract data
tt      = Psamp1.data(:,1)';
tmax    = max(tt);
ind     = 100e-3/(tt(2)-tt(1));          % Extract a part of the signal
 
tt      = tt(1:ind);
p1(1,:) = Psamp1.data(1:ind,2)'; % Mic 1 in coulmn 2 for rec. 1 - to p1
p1(2,:) = Psamp1.data(1:ind,3)';
p2(1,:) = Psamp2.data(1:ind,3)'; % Mic 1 in coulmn 2 for rec. 2 - to p1
p2(2,:) = Psamp2.data(1:ind,2)';                         
% Important to have control of these
 
% Plot for quality control and prep question
figure(11)
subplot(2,2,1)
plot1=plot(tt,p1)
plot2=plot(tt,p2)
xlabel('Time '), ylabel('Magnitude'), title( 'sumting ');
%xlim([0, 0.008]),ylim([-3,2]), legend([p1(1,1)], 'p1', 'p2')
set(gca,'fontsize',12,'fontweight','bold');
set(gcf,'units','centimeters','position',[2,1,29.7,21.0]) % Set size of plot to A4 size
 
%% Start computation 
% FFT
n  = 2^nextpow2( size(p1,2) );                          % Number of element in signal into fft 
fs = 44100;                              % (should be 2 to some power for better execution of fft)        
df = fs / n;

ff = fs*(0:(n-1))/n;            % Frequency vector
ww = 2*pi*ff;                   % Angular frequency

frecvec1 = fft(p1(1,:) ,n);     % Repeat for all data
frecvec2 = fft(p1(2,:) ,n);
frecvec3 = fft(p2(2,:) ,n);     % Repeat for all data
frecvec4 = fft(p2(1,:) ,n);

% Plot fft
figure(20) 
%frecvec = p1(1,:)
semilogx(ff,20*log10(frecvec1));
hold on;
semilogx(ff,20*log10(frecvec2));
semilogx(ff,20*log10(frecvec3));
semilogx(ff,20*log10(frecvec4));
title("FFT of P1");
xlim([f_l f_u]);
hold off;

%% Transfer function - see ISO standard for equations

H12I  = frecvec2 ./ frecvec1;           
H12II = frecvec3 ./ frecvec4;      

HC    = sqrt(H12I ./ H12II);    

H12   = H12I ./ HC;

figure(21)
semilogx(ff, abs(H12))
title("H12")

HI    = exp(-1i*ww/c  *s );
HR    = exp( 1i*ww/c  *s );          
 
%% Compute the absorption coefficient (alpha) and impedance (Z)
R =( (H12 - HI)./(HR - H12) ).*exp(2*1i*(ww/c)*x_1);    % eq. 17;

alpha = 1 - abs(R).^2;


Z = (1+R)./(1-R);
 
%% Plotting final results
 
figure(11)
subplot(2,2,2)
semilogx(ff,abs(R));
xlabel("Frequency [Hz]");
ylabel("Magnitude |R|");
title("Reflection Coefficient [|R|]")
%ylim
xlim([f_l f_u]);
 
figure(11)
subplot(2,2,3)
semilogx(ff,alpha);
grid on
xlim([f_l f_u]);
xlabel('Frequency (Hz)'), ylabel('Absorption coefficeint')
title('Absorption coefficient ')
 
figure(11)
subplot(2,2,4)
semilogx(ff,real(Z),'r')
hold on
semilogx(ff,imag(Z),'b');
semilogx(ff,abs(Z),'k');
grid on
% ylim([-4 4])
xlim([f_l f_u]);
title('Impedance Z')
xlabel('Frequency [Hz]' )
legend('')
hold off
 
 
%% Mikis model

%sigma = ; %Pa*s/m2
%e = ;
%f = ;
%w = ;
%k = ;
 
%Zc = ;
%r  = ;
%abs= ;
 
% Plot it
%figure
%semilogx(ff,aks)
%legend()
%title()
%xlabel()
%ylabel()
%xlim([f_l f_u]);
 
 
 
 %% Save the results
%  f = gcf;
%  %path2save = "C:\Users\erlen\Downloads\Lab2_TTT4250";
%  exportgraphics(f, ['ImpedanceTube.png'],'Resolution',450)
%  
%  
%  1+1   % Ez math?
 
 
 
 
 


