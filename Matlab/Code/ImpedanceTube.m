
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
 
%% Working frequency range 4.2
f_l =  100;          
f_u = .58*c/d;    
f_us = .45*c/s;
 
%% Load .etx files to Matlab
path2files = "C:\Users\erlen\Downloads\Lab2_TTT4250"
fileName1  = "C:\TTT4250_Lab2\Matlab\etx files\Imp_tube_12.etx"
fileName2  = "C:\TTT4250_Lab2\Matlab\etx files\Imp_tube_21.etx"

Tube_f1    = fileName1;
Tube_f2    = fileName2;
Psamp1     = importdata(Tube_f1,'\t',22); % 22 is the number of header lines
Psamp2     = importdata(Tube_f2,'\t',22);

%% Extract data
tt      = Psamp1.data(:,1)';
ind     = 150e-3/(tt(2)-tt(1));          % Extract a part of the signal
 
tt      = tt(1:ind);
p1(1,:) = Psamp1.data(1:ind,2)'; % Mic 1 in coulmn 2 for rec. 1 - to p1
p1(2,:) = Psamp1.data(1:ind,3)';
p2(1,:) = Psamp2.data(1:ind,3)'; % Mic 1 in coulmn 2 for rec. 2 - to p1
p2(2,:) = Psamp2.data(1:ind,2)';                         
% Important to have control of these
 
%% Plot for quality control and prep question
figure(11)
subplot(2,2,1)
plot1=plot(tt,p1)
hold on
plot2=plot(tt,p2)
ylim([-3 1.5])
legend("p1^I", "p2^I", "p2^{II}","p1^{II}" ,"Location", "best")
hold off
grid on
xlabel('Time [s] '), ylabel('Magnitude'), title( 'Impulse response ');
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

%% Plot fft
figure(20) 

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

HI    = exp(-1i*ww/c  *s );
HR    = exp( 1i*ww/c  *s );          
 
%% Compute the absorption coefficient (alpha) and impedance (Z)
R =( (H12 - HI)./(HR - H12) ).*exp(2*1i*(ww/c)*x_1);    % eq. 17;

alpha = 1 - abs(R).^2;
Z = (1+R) ./ (1-R) * rho * c;
 
%% Plotting final results
 
figure(11)
subplot(2,2,2)
semilogx(ff,20*log10(frecvec1));
hold on;
semilogx(ff,20*log10(frecvec2));
semilogx(ff,20*log10(frecvec3));
semilogx(ff,20*log10(frecvec4));
xlabel("Frequency [Hz]");
ylabel("Magnitude");
title("FFT of the Impulse respons")
set(gca,'fontsize',12,'fontweight','bold');
legend("FFT(p1^I)","FFT(p2^I)","FFT(p2^{II})" ,"FFT(p1^{II})","Location", "best")
grid on
%ylim
xlim([f_l f_u]);
hold off
 
figure(11)
subplot(2,2,3)
semilogx(ff,alpha);
hold on
semilogx(ff,abs(R));
grid on
legend("Absorption", "Reflection", "Location", "best")
hold off
xlim([f_l f_u]);
xlabel('Frequency [Hz]'), ylabel('Magnitude')
title('Absorption-and Reflection coefficient ')
set(gca,'fontsize',12,'fontweight','bold');


figure(11)
subplot(2,2,4)
semilogx(ff,real(Z),'r')
hold on
semilogx(ff,imag(Z),'b');
semilogx(ff,abs(Z),'k');
grid on
% ylim([-4 4])
xlim([f_l f_u]);
title('Acoustic Impedance Z')
xlabel('Frequency [Hz]' )
legend('Re[Z]', "Im[Z]", "|Z|", "Location", "best")
set(gca,'fontsize',12,'fontweight','bold');
hold off
 
%% Mikis model

f = [100:1:2000];
omega = 2*pi*f;

rho_0 = 1.225;      % [Kg.m-3] density at rest of air at 18C, 1atm
c_0   = 342.2;      % [m.s-1] speed of sound in air at 18C, 1atm
P_0   = 1.0132e+05; % [N.m-2] atmospheric pressure at 18C, 1atm
sigma = 9100       % [N.s.m-4]?static air flow resistivity of material
h     = 0.1        % [m]?thickness of material
X = f/sigma;

Z_DB70_Mik90 = rho_0*c_0*( 1 + 5.50*(X*1000).^(-0.632) ...
                            - i*8.43*(X*1000).^(-0.632) ); 

k_DB70_Mik90 = omega/c_0 .* (-i) .* ( 11.41*(X*1000).^(-0.618) ...
                                      + i* (1 + 7.81*(X*1000).^(-0.618) ) );

Z = Z_DB70_Mik90;
Z = -1i.*Z_DB70_Mik90./tan(k_DB70_Mik90*h);
Z_cot = -1i.*Z_DB70_Mik90./cot(k_DB70_Mik90*h);

R = (Z-rho_0*c_0)./(Z+rho_0*c_0);
R_cot = (Z_cot-rho_0*c_0)./(Z_cot+rho_0*c_0);
a = 1-abs(R).^2
a_cot = 1-abs(R_cot).^2
figure(10)
subplot(1,2,2)
%semilogx(f, real(Z))
semilogx(f, abs(Z_cot))
hold on
%semilogx(f, imag(Z))
semilogx(f, abs(Z))
grid on
title('Specific Impedance Z_C')
xlabel('Frequency [Hz]' )
ylabel("Magnitude")
%legend('Re[Z]', "Im[Z]", "|Z|")
legend("tan", "cot")
set(gca,'fontsize',12,'fontweight','bold');
set(gcf,'units','centimeters','position',[2,1,29.7,11.0])
hold off

figure(10)
subplot(1,2,1)
semilogx(f, abs(a_cot))
hold on
semilogx(f, abs(a))
title('Absorption-and Reflection Coefficient')
xlabel('Frequency [Hz]')
ylabel("Magnitude")
%legend('Absorption', "Reflection", "Location", "best")
legend("cot", "tan")
set(gca,'fontsize',12,'fontweight','bold');
grid on
hold off

%xlim([f_l f_u]);
%% Save the results
exportgraphics(figure(10), ['Mikis_Model.png'],'Resolution',450)
exportgraphics(figure(11), ['Impedance_Tube.png'],'Resolution',450)
 

%  1+1   % Ez math?
 
 
 
 
 


