function hw8p3
e = 1.60217662e-19; % Electron Charge [J]
r_inner = .10/2;    % Inner Radius [m]
r_outer = .15/2;    % Outer Radius [m]
Ae = (pi*r_outer^2)-(pi*r_inner^2); % Exit Area [m^2]
ni = 5e17;      % Ion Plasma Density [m^-3]
B = 200;        % Radial Magnetic Field [G]
m =9.10938e-31; % Electron Mass [kg]
M = 2.18e-25;   % Ion Mass [kg]
Te = 20;        % Electron Temperature [eV]
Vd = 300;       % Discharge Champer Potential [V]

k = 1.38064852e-23; % Boltzmann Constant [J/K]
k2 = 8.6173303e-5;  % Boltzmann Constant [eV/K]

%Part A
Ii = ni*e * sqrt((2*e*Vd)/(M))*Ae;  % Beam Current [A]
P = Ii*Vd;                          % Beam Power [W]

%Part B
r_L = ((m*1000)/(e*B))*sqrt((8*(k/k2)*Te)/(pi*m))*1000*10; % Larmor Radius [mm]

%Part C
wB = (e*B)/(m);                     % Cyclotron Frequency [1/s]
Q = 6.5e-13/((3/2)*Te)^2;           % Collisional Cross Section
vth = sqrt((8*(k/k2)*Te)/(pi*m));   % Collision Speed
nu = ni*Q*vth;                      % Collision Frequency [1/s]

omega = sqrt(wB^2/nu^2); % Hall Parameter

%Part D
gamma= 0.9;     % Thrust Coefficient Factor
eta_m = 0.8;    % Mass Utilization Efficiency


mdot_i = (Ii*M)/e;
vi = sqrt((2*e*Vd)/(M));
Isp = (gamma*eta_m*vi)/9.81;
T = gamma*mdot_i*vi;

%Part E
IH = ni*e*((.15-.1)/2)*(Vd/(B*10^-4));
I_H = (Ii/(2*pi*((r_inner+r_outer)/2)*(B*10^-4)))*sqrt((M*Vd)/(2*e));
end
