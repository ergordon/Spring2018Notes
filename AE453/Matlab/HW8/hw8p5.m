function hw8p5
clc; clear;
m =9.10938e-31;
M = 2.18e-25;           % Ion Mass [kg]
k = 1.3806e-23;         % Boltzmann Constant [J/K]
k2 = 8.6173e-5;         % Boltzmann Constant [eV/K]
K = (k/k2);             % Boltzmann Constant [J/eV]
n_o = 10^13 * (1/1e-6); % Neutral Xenon Gas Density [m^-3]
Volume = 10^4*1e-6;     % Plasma Volume [m^3]
Area = 200*0.0001;      % Ion Loss Area [m^2]

EE = (2*n_o*Volume)/Area;

eo = 8.854e-12;         % Vacuum Dielectric Constant [C^2/Nm^2]
h = 6.6262e-34;         % Planks Constant [Js]
qe = 1.60217662e-19;    % Electron Charge [C]
ao = (eo*h^2)/(pi*m*qe^2);  % Atomic Cross Section [m^2]
Q = pi*ao^2*4.38;           % Emperical Cross Section [m^2]

%% This approach uses equations fitted to the empiracle ionization distributions to calculate the maxwellian distribution given a specific value of Te.

Te_low = 1;
Te_high = 3;
while Te_high-Te_low > .0000000000001
    Te = (Te_low+Te_high)/2;
    if Te < 5
        %Maxwellian Ionization Rate Coefficient
        MOVi = 10^-20 * ((3.97+0.643*Te - 0.0368*Te^2)*exp(-12.127/Te))*sqrt((8*K*Te)/(pi*m));
    else
        %Maxwellian Ionization Rate Coefficient
        MOVi = 10^-20 * (-(1.031e-4*Te^2) + 6.386*exp(-12.127/Te))*sqrt((8*K*Te)/(pi*m));
    end

    Vth = sqrt((8*K*Te)/(pi*m));
    OVi = 0.0005*Q*Vth + 0.9995*MOVi;
    Te_new = (EE^2 * OVi^2)*(M/K);
    
    if Te_new >= Te
        Te_high = Te;
    else
        Te_low = Te;
    end
end

%% This approach assumes a linear relation between the maxwellian ionizaion reaction rates and solves for a new ionization coefficient rate. With this,Te is recalculated and compared to the guessed Te value.

Te_low = 1.5;
Te_high = 2;
OVi_low = 1.16e-17;
OVi_high = 1.08e-16;

while Te_high-Te_low > .000000000001
    ovi = (OVi_high+OVi_low)/2;
    t = (Te_high+Te_low)/2;
    
    Vth = sqrt((8*K*Te)/(pi*m));
    OVI = 0.0005*Q*Vth + 0.9995*ovi;
    Te_new = (EE^2 * OVI^2)*(M/K);
    
    if Te_new >= t
        OVi_high = ovi;
        Te_high = t;
    else
        OVi_low = ovi;
        Te_low = t;
    end
end

end