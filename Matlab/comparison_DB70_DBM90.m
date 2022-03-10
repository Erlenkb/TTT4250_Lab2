% comparison of models by Delany-Bazley 
% and Delany-Bazley-Miki
%
% M. E. Delany and E. N. Bazley, 
% Acoustical properties of fibrous absorbent materials,
% Applied Acoustics (3), 1970, pp. 105-116
%
% Y. Miki
% Acoustical properties of porous materials 
% - modifications of Delany-Bazley models -
% J. Acoust. Soc. Jpn. (E), 11 (1), 1990, pp. 19-24
%
% Copyleft 2006 luc.jaouen@matelys.com
% cf. APMR on the web,
% PropagationModels/MotionlessSkeleton/DelanyBazleyMikiModel.html
% for more information

close all
clear all

f = [100:1:2000];
omega = 2*pi*f;

rho_0 = 1.225;      % [Kg.m-3] density at rest of air at 18C, 1atm
c_0   = 342.2;      % [m.s-1] speed of sound in air at 18C, 1atm
P_0   = 1.0132e+05; % [N.m-2] atmospheric pressure at 18C, 1atm
sigma = 9100       % [N.s.m-4] static air flow resistivity of material
h     = 0.1        % [m] thickness of material


%%%%%
%%%%% Compute variable X and print frequency
%%%%% limits of validity for the two models 
%%%%%
X = f/sigma;
f_min = 0.01*sigma
f_max = 1.00*sigma


%%%%%
%%%%% Delany and Bazley model
%%%%% (NB: gamma = alpha + j beta = j k )
%%%%%
Z_DB70 = rho_0*c_0*( 1 + 9.08*(X*1000).^(-0.75) ...
                     - i*11.9*(X*1000).^(-0.73) ); 

k_DB70 = omega/c_0 .* (-i) .* ( 10.3*(X*1000).^(-0.59) ...
                                + i* ( 1 + 10.8*(X*1000).^(-0.70) ) );

K_DB70 = Z_DB70.*omega./k_DB70;
rho_DB70 = k_DB70.*Z_DB70./omega;


%%%%%
%%%%% Revised expressions of Delany and Bazley model by Miki 
%%%%% (NB: gamma = alpha + j beta = j k )
%%%%%

Z_DB70_Mik90 = rho_0*c_0*( 1 + 5.50*(X*1000).^(-0.632) ...
                            - i*8.43*(X*1000).^(-0.632) ); 

k_DB70_Mik90 = omega/c_0 .* (-i) .* ( 11.41*(X*1000).^(-0.618) ...
                                      + i* (1 + 7.81*(X*1000).^(-0.618) ) );

Z = Z_DB70_Mik90;
Z = -1i.*Z_DB70_Mik90./tan(k_DB70_Mik90*h);

R = (Z-rho_0*c_0)./(Z+rho_0*c_0);
a = 1-abs(R).^2

figure(10)
subplot(1,2,2)
semilogx(f, real(Z))
hold on
semilogx(f, imag(Z))
semilogx(f, abs(Z))
grid on
title('Impedance Z')
xlabel('Frequency [Hz]' )
ylabel("Magnitude")
legend('Re[Z]', "Im[Z]", "|Z|")
set(gca,'fontsize',12,'fontweight','bold');
set(gcf,'units','centimeters','position',[2,1,29.7,11.0])
hold off

figure(10)
subplot(1,2,1)
semilogx(f, abs(R))
hold on
semilogx(f, abs(a))
title('Absorption and Reflection')
xlabel('Frequency [Hz]')
ylabel("Magnitude")
legend('Absorption coefficient', "Reflection coefficient", "Location", "best")
set(gca,'fontsize',12,'fontweight','bold');
grid on
hold off

K_DB70_Mik90 = Z_DB70_Mik90.*omega./k_DB70_Mik90;

figure(32)
semilogx(f,K_DB70_Mik90)

rho_DB70_Mik90 = k_DB70_Mik90.*Z_DB70_Mik90./omega;


%%%%%
%%%%% Compute sound absorption using the two models 
%%%%% for a sample of thickness d backed by a rigid 
%%%%% and impervious wall under at room temperature
%%%%% and pressure conditions 
%%%%%

Z = -j.*Z_DB70./tan(k_DB70*h);
alpha_DB70 = 1 - ( abs( (Z-rho_0*c_0)./(Z+rho_0*c_0) ) ).^2;


Z = -j.*Z_DB70_Mik90./tan(k_DB70_Mik90*h);
alpha_DB70_Mik90 = 1 - ( abs( (Z-rho_0*c_0)./(Z+rho_0*c_0) ) ).^2;


%%%%%
%%%%% Compare results
%%%%%

figure(1)
set(gca,'FontSize',16)
plot(f,real(rho_DB70)/rho_0,'k-','LineWidth',2)
hold on
plot(f,imag(rho_DB70)/rho_0,'b-','LineWidth',2)
plot(f,real(rho_DB70_Mik90)/rho_0,'r--','LineWidth',2)
plot(f,imag(rho_DB70_Mik90)/rho_0,'m--','LineWidth',2)
xlabel('Frequency (Hz)')
ylabel('Normalized dynamic density')
legend('DB70 (Re)','DB70 (Im)','DB70+Mik90 (Re)','DB70+Mik90 (Im)')

figure(2)
set(gca,'FontSize',16)
plot(f,real(K_DB70)/P_0,'k-','LineWidth',2);
hold on
plot(f,imag(K_DB70)/P_0,'b-','LineWidth',2);
plot(f,real(K_DB70_Mik90)/P_0,'r--','LineWidth',2);
plot(f,imag(K_DB70_Mik90)/P_0,'m--','LineWidth',2);
xlabel('Frequency (Hz)')
ylabel('Normalized dynamic bulk modulus')
legend('DB70 (Re)','DB70 (Im)','DB70+Mik90 (Re)','DB70+Mik90 (Im)')

figure(3)
set(gca,'FontSize',16)
semilogx(f,alpha_DB70,'k-','LineWidth',2);
hold on
semilogx(f,alpha_DB70_Mik90,'r--','LineWidth',2);
xlim([100 2000])
grid on
xlabel('Frequency (Hz)')
ylabel('Sound absorption coefficient')
legend('DB70','DB70+Mik90')
