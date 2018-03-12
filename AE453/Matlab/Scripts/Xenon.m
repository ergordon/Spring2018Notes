function Xenon
%Xenon Xe
MW = 131.293; %Molecular Weight [amu]
Qxe1 = 1.56e-19; %Charge of Xenon Ion [c]
Mxe1 = 2.18e-25; %Mass of Xenon (Ion) [kg]
epsilon_i = 12.13; %Ionization Energy wrt. Ground State [eV]

specificHeat(2050)
end

function Cp = specificHeat(T)
if (T>298)&&(T<6001)
    T=T/1000;
    A = 20.78600;
    B = 7.449320e-7;
    C = -2.049401e-7;
    D = 1.066661e-8;	
    E = 2.500261e-8;
else
    Cp = 0;    
end
Cp = A+(B*T)+(C*(T^2))+(D*(T^3))+(E/(T^2));
end