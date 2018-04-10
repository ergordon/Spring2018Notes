clear; clc; format short
%% Define Parameters

%Given
mdot = 61; % Mass Flow Rate [mg/s]
AeAt = 33; % Nozzle Expansion Ratio
T = 2050;  % Peak Temperature [K]
Pc = 200;   % Feed Pressure [psia]
Pin = 500; % Input Power [kg/s]
x = 0.7;   % Dissasociation
gamma = 1.351; % Adiabatic Index (See Notes)

%Universal
MW_N2H4 = 32.0452; % Molecular Weight Hydrazine [g/mol]
MW_NH3 = 17.03052; % Molecular Weight Ammonia [g/mol]
MW_N2 = 28.014;   % Molecular Weight N2 [g/mol]
MW_H2 = 2.016;   % Molecular Weight H2 [g/mol]

R = 8.314; % Universal Gas Constant [J/mol.K]

%% Conversions

mdot = mdot*1e-6; % [mg/s] to [kg/s]
Pc = Pc*6894.76;  % [psia] to [Pa]=[kg/m.s^2]

%% Find Weighted Averages

mBar_NH3 = (4/3)*(1-x); % Moles of NH3 [mol]
mBar_N2 = (1/3)*(1+2*x); % Moles of N2 [mol]
mBar_H2 = 2*x; % Moles of H2 [mol]
mBar_all = mBar_NH3+mBar_N2+mBar_H2;

mFrac_NH3 = mBar_NH3/mBar_all; % Mass of NH3 [g]
mFrac_N2 = mBar_N2/mBar_all; % Mass of N2 [g]
mFrac_H2 = mBar_H2/mBar_all; % Mass of H2 [g]

%% Weighted Average for Molecular Weight [g/mol]
MW_avg = (mFrac_NH3*MW_NH3)+(mFrac_N2*MW_N2)+(mFrac_H2*MW_H2);

%% Specific Gas Constant [J/g.K]
R_N2H4 = R/MW_avg;    % [J/g.K]
R_N2H4 = R_N2H4*1000; % [J/kg.K]

%% Equation 7.59
syms Pe
% Exit Area Pressure [Pa]
Pe = solve((1/AeAt) == ((gamma+1)/(2))^((1)/(gamma-1))*(Pe/Pc)^(1/gamma)*sqrt(((gamma+1)/(gamma-1))*(1-(Pe/Pc)^((gamma-1)/gamma))),Pe);

%% Equation 7.58
% Throat Area [m^2]
A_t = mdot/(Pc*sqrt((gamma/(R_N2H4*T))*(2/(gamma+1))^((gamma+1)/(gamma-1))));
% Throat Area [mm^2]
A_t = A_t*1000^2;
% Exit Area [mm^2]
A_e = A_t*AeAt;

%% Equation 7.48
% Exit Velocity [m/s]
Ue = sqrt(((2*gamma*R_N2H4*T)/(gamma-1))*(1-((Pe)/(Pc))^((gamma-1)/(gamma))));

% Thrust [N]
Thrust = mdot*Ue+A_e*10e-6*Pe;

% Specific Impulse [sec]
I_SP = Ue/9.81;

% Thruster Efficiency [%]
eta = 0.5*Thrust*Ue / 500;

sprintf('The Aerojet MR-501B has: \n A Thrust of %d (mN);\n Specific impulse %d (sec);\n A Thrust Efficiency of %0.2f (percent); \n Throat Diameter %0.5f (mm);\n Nozzle Exit Diameters %0.5f (mm)', Thrust*1000,I_SP, eta*100, vpa(sqrt(A_t/pi)), vpa(sqrt(A_e/pi)))

%No Resisitve Heating
Cp = gamma*R/(gamma-1)
Ue = sqrt(Cp*2*T)
T = mdot*Ue*1000
ISP = Ue/9.81