function hw7a
clc; clear;

MW = 39.948;        % Argon Molecular Weight [g/mol]
l = 0.05;           % Channel Length [m]
area = 0.03*0.03;   % Area of Channel
volume = 0.05*area; % Volume of Channel
R_A = 208.13;       % Specific Gas Constant for Argon [J/kgK]
x = [0:0.001:l];    % x values along the channel
P_o = 66.661185;    % Peak Pressure [Pa]
T_o = 273;          % Peak Temperature [K]
rho_o = P_o / (R_A * T_o) % Total Density of Heavy Particles
J = 20e3;                 % Initial Current Sheet Disharge [A]
Lprime = 0.6E-6;          % Distributed Inductance [H]

%Part 1: Density Profiles
for i=1:length(x)
    rhoCD(i) = 1*rho_o;                   % Constant Density
	rhoLDD(i) = (1-(x(i)/l))*rho_o;       % Linear Decreasing Density
	rhoSQD(i) = sqrt((1-(x(i)/l)))*rho_o; % Square Root Dependence, Decreasing Density
end

%Part 2: Total Mass Injected [kg] for
syms xx(t)
tmCD = area*rho_o*l                                   % Constant Density
tmLDD = area*rho_o*(l-((l^2)/(2*l)))                  % Linear Decreasing Density
tmSQD = area*rho_o*integral(@(xx) sqrt(1-(xx/l)),0,l) % Square Root Dependence, Decreasing Density

 
%Part 1 Solution
figure(1)
plot(x,rhoCD,'linewidth',2)
hold on
plot(x, rhoLDD,'linewidth',2)
hold on
plot(x, rhoSQD,'linewidth',2)
    grid on
    grid minor
    xlabel 'Location Along Channel Length [m]'
    ylabel 'Density [Torr]'
    title 'Particle Density'
    legend('Constant Density','Linear Decreasing Density','Square Root Dependence,Decreasing Density','location','bestoutside')

    
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartA.pdf')

%Part 2 Plot
x = [0:0.001:l];
for i=1:length(x)
    mCD(i) = area*rho_o*x(i) % Constant Density
    mLDD(i) = area*rho_o*(x(i)-((x(i)^2)/(2*l))) % Linear Decreasing Density
    xx = x(i)
    mSQD(i) = area*rho_o*integral(@(xx) sqrt(1-(xx/l)),0,xx) % Square Root Dependence, Decreasing Density
end
 
figure(2)
plot(x,mCD,'linewidth',2)
hold on
plot(x,mLDD,'linewidth',2)
hold on
plot(x,mSQD,'linewidth',2)
    grid on
    grid minor
    xlabel 'Location Along Channel Length [m]'
    ylabel 'Mass Accumulated [kg]'
    title 'Total Mass Injected'
    legend('Constant Density','Linear Decreasing Density','Square Root Dependence,Decreasing Density','location','bestoutside')

    
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartB.pdf')
end