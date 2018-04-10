clear; clc; format short
%% Define Parameters

%Given
Pin = 100; % Input Power [kW]
Tc = 4000; % Chamber Temperature [K]
Pc = 1; % Chamber Pressure [atm]
mdot = 220; % Mass Flow Rate [mg/s]
Pe = 0.01; % Exit Pressure [atm]

%Universal
m_H2 = 1.673723e-27; %MW of H2 [kg]
R = 8.314; % Universal Gas Constant [J/mol.K]

%% Conversions

mdot = mdot*1e-6;  % [mg/s] to [kg/s]
Pc = Pc*101325;    % [atm] to [Pa]=[kg/m.s^2]
Pe = Pe*101325;    % [atm] to [Pa]=[kg/m.s^2]
Pin = Pin*1000;    % [kW] to [W]

%% Dissassociation Calculation
% At Tc = 4000 [K]
K = exp(0.934);

X2HXH2 = K/(Pc/Pe);

sol = roots([-1, -X2HXH2, X2HXH2]);

X = sol(2); % Dissasociation
sprintf('Equillibrium Degree of Dissasociatio in the Combustion Chamber is %0.3f(percent)', X*100)

No = 1/m_H2;
alpha2 = 1-X;
alpha1 = 2*X;

%% Combustion Chamber Enthalpy and Ratio of Specific Heat
hc = No*(alpha2*((9/2)*K*Tc)+alpha1*((5/2)*K*Tc+0.5*X));
sprintf('The Combustion Chamber Enthalpy is %d', hc)

syms gamma
eqn1 = (hc/Tc) == (gamma*R)/((gamma-1)*m_H2);

gamma = vpa(solve(eqn1,gamma));
sprintf('The Ratio of Specific Heats is %0.3f', gamma)

%% Area Calculations
% Equation 7.58 - Throat Area [m^2]
A_t = mdot/(Pc*sqrt((gamma/(R*Tc))*(2/(gamma+1))^((gamma+1)/(gamma-1))));
% Equation 7.59 - Exit Area [m^2]
syms A_e
A_e = solve((A_t/A_e) == ((gamma+1)/(2))^((1)/(gamma-1))*(Pe/Pc)^(1/gamma)*sqrt(((gamma+1)/(gamma-1))*(1-(Pe/Pc)^((gamma-1)/gamma))),A_e);

% Throat Area [mm^2]
A_t = A_t*1000^2;
% Exit Area [mm^2]
A_e = A_e*1000^2;

% Throat Diameter [mm]
D_t = vpa(sqrt(A_t/pi));
% Exit Diameter [mm]
D_e = vpa(sqrt(A_e/pi));

%% Part C Answers
% Equation 7.48 - Exit Velocity [m/s]
Ue = sqrt(((2*gamma*R*Tc)/(gamma-1))*(1-((Pe)/(Pc))^((gamma-1)/(gamma))));
% Thrust [N]
Thrust = mdot*Ue+A_e*10e-6*Pe;
% Specific Impulse [sec]
I_SP = Ue/9.81;
% Thruster Efficiency [%]
eta = 0.5*Thrust*Ue / Pin;

sprintf('Assuming Frozen Flow with this Chamber Composition we get:\n A Thrust of %0.3f (mN);\n Specific impulse %0.2f (sec);\n A Thrust Efficiency of %0.2f (percent); \n Throat Diameter %0.5f (mm);\n Nozzle Exit Diameters %0.5f (mm)', Thrust*1000,I_SP, eta*100, D_t, D_e)

%% Part D Answers

% Isentropic Pressure Relation - T2/T1 = (P2/P1)^(1-(1/gamma))

% Exit Temperature [K]
Te = Tc*(Pe/Pc)^((gamma-1)/gamma);
sprintf('Exit Temperature for an Isentropic Nozzle is %d [K]',Te)
% Dissassociation Calculation - At Te = 150 [K]
K = exp(-164.005);

X2HXH2 = K/(Pc/Pe);

sol = roots([-1, -X2HXH2, X2HXH2]);

X = sol(1); % Dissasociation
sprintf('Equillibrium Degree of Dissasociation in the Combustion Chamber is %0.3f(percent)', X*100)

No = 1/m_H2;
alpha2 = 1-X;
alpha1 = 2*X;
xi_Frozen = 1;
xi_Equillibrium = 0;

hc = No*(alpha2*((9/2)*K*Tc)+alpha1*((5/2)*K*Tc+0.5*X));
he = 0.5*xi_Equillibrium*alpha1*No*X
sprintf('The Exit Enthalpy is %0.5f', he)

% Equation 7.48 - Exit Velocity [m/s]
%Ue = sqrt(((2*gamma*R*Te)/(gamma-1))*(1-((Pe)/(Pc))^((gamma-1)/(gamma))))
%1D Energy Equation for Nozzle Flow
%Ue = sqrt(2*(hc-he))
%Equation 7.61
Ut = sqrt(gamma*R*Tc);
Ue = Ut*sqrt(((gamma+1)/(gamma-1))*(1-((Pe)/(Pc))^((gamma-1)/(gamma))));

% Thrust [N]
Thrust = mdot*Ue+A_e*10e-6*Pe;
% Specific Impulse [sec]
I_SP = Ue/9.81;
% Thruster Efficiency [%]
eta = 0.5*Thrust*Ue / Pin;

sprintf('Assuming The Flow Maintains Equillibrium We Get:\n A Thrust of %0.3f (mN);\n Specific impulse %0.2f (sec);\n A Thrust Efficiency of %0.2f (percent)', Thrust*1000,I_SP, eta*100)
