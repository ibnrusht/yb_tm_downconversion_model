clear;
%constants
c = 3*10^10; % сm/sec - light velocity
Na = 6.02*10^23; % 1/mol Avogadro constant
p = 6.09; % g/cm3 P. Lecoq M. Schussler "progress and prospect in the development of new scintillators for future"
MM = 256; % g/mol - molar mass
l = 688*10^-7; % cm - pump wavelength
Spp = 1.5*10^-21; % сm^2 - absorption cross section было 1.5*10^-21
Sle = 4.2*10^-21; % сm^2 - laser emission cross-section on 1050 nm (было 4.2*10^-21)
Sla = 0.057*10^-21; % cm^2 - laser absorption cross-section on 1050nm (2009 / Vol. 17, No. 20 / OPTICS EXPRESS, Direct Comparison of Yb3+)
h_plank = 6.63*10^-34; % J/sec - Planck constant
k_b = 1.38*10^(-23); % J/K - Boltzmann constant
k_Yb = 0.2;
Nyb = k_Yb*(p*Na)/MM; % 1/сm3 - Yb concentration
Ntm = 0.002*(p*Na)/MM; % 1/сm3 - Tm concentration
ti = 1e-8; % s - puls duration 0.01
T = 0.1; % s - pulse period 0.06
k = 0.4; % Tm ratio
A_Tm_Yb = 2.5*10^-20;
A_Tm_Yb_2 = 6.5*10^-19;
ph = 428; % phonon energy cm-1


N0=[k*Ntm;0;0;Nyb;0;(1-k)*Ntm;0;0]; % initial conditions
% N0=[Ntm;0;0;Nyb;0];

time=0:1e-7:0.015;

start = 80;
W = [];

for i=0:20:20
    temp = start + i;
    [t,n]=ode23(@diff_eq,time,N0,[],temp,A_Tm_Yb,A_Tm_Yb_2,ph);
    writematrix(n(:,2) + n(:,7), strcat(num2str(k_Yb),'Yb_', num2str(temp), 'K_3H4.txt'))
    % writematrix(n(:,2), strcat(num2str(k_Yb),'Yb_', num2str(temp), 'K_3H4.txt'))
    dE47 = 2800;
    C47 = (1-exp(-ph*h_plank*c/(k_b*temp)))^round(-dE47/ph);
    W = [W, 1/(A_Tm_Yb_2*C47*Nyb)];
end

writematrix(W, strcat(num2str(k_Yb),'Yb_K_3H4.txt'))
delta = n(end,1) + n(end,6) - Ntm
% delta = n(end, 1) - Ntm
