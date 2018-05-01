function hw8p1
r = 0.05; %Chamber Radius [m]
l = 0.15; % Chamber Length [m]
A = pi*r^2; %Grid Area [m^2]
A_a = pi*r^2 + 2*pi*r*l; %Anode Area [m^2]
T_g = 0.8; %Grid Transparency
n_o = 10^18; %Neutral Particle Density [m^-3]
V = A*l; % Volume [m^3]
Ui = 12.13;% Ionization Potential [eV]
Ue = 10; % Excitation Potential [eV]
m =9.10938e-31;
M = 2.18e-25;

k = 1.38064852e-23; % Boltzmann Constant [J/K]
k2 = 8.6173303e-5; % Boltzmann Constant [eV/K]

Te = [3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]; % Electron Temperature [eV]
OVi = [1.08e-15 2.13e-15 3.59e-15 5.43e-15 7.61e-15 1.01e-14 1.28e-14 1.57e-14 1.88e-14 2.20e-14 2.53e-14 2.86e-14 3.20e-14 3.55e-14 3.90e-14];% Ionization Rate Coefficient [m^3/s]
OVe = [2.66e-15 4.66e-15 7.12e-15 9.93e-15 1.30e-14 1.61e-14 1.94e-14 2.26e-14 2.57e-14 2.87e-14 3.14e-14 3.34e-14 3.41e-14 3.21e-14 2.48e-14];% Excitation Rate Coefficient [m^3/s]

for i = 1:length(Te)
    v_a = sqrt(((k/k2)*Te(i))/M);
    eta_d(i) = (((2*n_o*OVi(i)*V)/(T_g*A*v_a))*(Ui + (OVe(i)/OVi(i))*Ue))+((1/(T_g))*(2.5*Te(i) + 2*Te(i)*log((A_a/A)*sqrt((2*M)/(m*pi)))));
end

plot(Te,eta_d,'linewidth',2)
title('Discharge Loss V. Electron Temperature')
ylabel('Discharge Loss (eV/ion)')
xlabel('Electron Temperature [eV]')
grid on
grid minor

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
%print(gcf,'-dpdf','Problem1.pdf')
end

