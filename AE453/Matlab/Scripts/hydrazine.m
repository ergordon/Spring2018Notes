function hydrazine
%Hyrdrazine H4N2
MW = 32.0452; %Molecular Weight [amu]

specificHeat(1000)
end

function Cp = specificHeat(T)
if (T>800)&&(T<2001)
    T=T/1000;
    A = 35.1824;
    B = 96.05260;
    C = -40.50130;
    D = 6.668070;	
    E = -0.874233;
elseif (T>2000)&&(T<6000)
    T=T/1000;
    A = 121.4010;
    B = 4.816880;
    C = -0.763012;
    D = 0.043232;
    E = -40.78650;
else
    Cp = 0;    
end
Cp = A+(B*T)+(C*(T^2))+(D*(T^3))+(E/(T^2));
end