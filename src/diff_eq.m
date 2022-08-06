function dndt=diff_eq(t,n,temp,A_Tm_Yb,A_Tm_Yb_2,ph) % îïèñàíèå ñèñòåìû êèíåòè÷åñêèõ óðàâíåíèé

h_plank = 6.63*10^-34;
c = 3*10^10;
S_abs = 2.*10^-20; % absorption cross-section of Tm transition 3H6-3F3
P = 0.0001;
lambda = 688*10^-7;
S = 2*10^-3;
tn = 5e-8;
ti = 5e-9;

k = 1.38*10^(-23); % J/K - Boltzmann constant
T = temp; % K - crystal temperature

A_3H4 = 400;
A_3F3 = 1343.99;
A_2F52 = 400;
A47 = A_Tm_Yb;
A472 = A_Tm_Yb_2;
MPR = 260000;
dE54 = 1700;
dE47 = 2300;
dE472 = 2500;
phonon = ph;

C54 = (1-exp(-phonon*h_plank*c/(k*T)))^round(-dE54/phonon);
C47 = (1-exp(-phonon*h_plank*c/(k*T)))^round(-dE47/phonon);
C472 = (1-exp(-phonon*h_plank*c/(k*T)))^round(-dE472/phonon);

dN_3H6dt = -P*exp(-((t-tn)^2/(2*ti^2))^4)*S_abs*lambda*n(1)/(h_plank*c*S) + A_3H4*n(2) + A_3F3*n(3) + A47*C472*n(2)*n(4);
dN_3H4dt = - A_3H4*n(2) + C54*MPR*n(3) - A47*C472*n(2)*n(4);
dN_3F3dt = P*exp(-((t-tn)^2/(2*ti^2))^4)*S_abs*lambda*n(1)/(h_plank*c*S) - (C54*MPR + A_3F3)*n(3);
dN_2F72dt = - A47*C472*n(2)*n(4) - A472*C47*n(7)*n(4) + A_2F52*n(5);
dN_2F52dt = A47*C472*n(2)*n(4) + A472*C47*n(7)*n(4) - A_2F52*n(5);

dN_3H62dt = -P*exp(-((t-tn)^2/(2*ti^2))^4)*S_abs*lambda*n(6)/(h_plank*c*S) + A_3H4*n(7) + A_3F3*n(8) + A472*C47*n(7)*n(4);
dN_3H42dt = - A_3H4*n(7) + C54*MPR*n(8) - A472*C47*n(7)*n(4);
dN_3F32dt = P*exp(-((t-tn)^2/(2*ti^2))^4)*S_abs*lambda*n(6)/(h_plank*c*S) - (C54*MPR + A_3F3)*n(8);

dndt=[dN_3H6dt;dN_3H4dt;dN_3F3dt;dN_2F72dt;dN_2F52dt;dN_3H62dt;dN_3H42dt;dN_3F32dt];
