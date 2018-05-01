function hw8p2
clc;clear;

M = 2.18e-25;   % Ion Mass [kg]
n_o = 10^13 * (1/1e-6); %Beutral Xenon Gass Density [m^-3]
Volume = 10^4*1e-6; % Plasma Volume [m^3]
Area = 200*0.0001; % Ion Loss Area [m^2]

k = 1.38064852e-23; % Boltzmann Constant [J/K]
k2 = 8.6173303e-5;  % Boltzmann Constant [eV/K]

K = (k/k2);
EE = (2*n_o*Volume)/Area;

OVi_low = 1.08e-16;
OVi_high = 1.2775e-16; % From Hand Calculations
Te_low = 2;
Te_high = 2.03125; % From Hand Calculations

%Assuming linear relation between ionization rates.
while Te_high-Te_low > .00000001
    ovi = (OVi_high+OVi_low)/2;
    t = (Te_high+Te_low)/2;
    Te_new = (EE^2 * ovi^2)*(M/K);
    if Te_new >= t
        OVi_high = ovi;
        Te_high = t;
    else
        OVi_low = ovi;
        Te_low = t;
    end
end

Te_guess = Te_low;
OVi = OVi_high
Te_calculated = (EE^2 * OVi^2)*(M/K)

end