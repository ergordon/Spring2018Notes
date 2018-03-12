function hw4p4
    %SAHA(1,5000)
    %%{
    T = [2000:100:6000];
    P = [10^-3, 10^-1, 10];
    for i=1:length(T)
        AlphaPT(1,i)=SAHA(P(1),T(i));
        AlphaPT(2,i)=SAHA(P(2),T(i));
        AlphaPT(3,i)=SAHA(P(3),T(i));
    end
    plot(T,AlphaPT,'LineWidth',2)
    title('Ionization Dependence of Xenon on Temperature and Pressure', 'FontSize', 18);
ylabel('Degree of Ionization \alpha', 'FontSize', 13)
xlabel('Temperature [K]', 'FontSize', 13)
grid on
grid minor
legend('P = 10^{-3} Torr',...
    'P = 10^{-1} Torr',...
    'P = 210 Torr',...
    'Location','southeast')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','p4.pdf');

    
    %}

    %Atomic Xenon Partition Function (AXPF)
    function f_A = AXPF(T)
        th = T/1; %Dimensionless Temperature T/To
        f_A = -1.08e-21*th^5 + 1.86e-16*th^4 -6.49e-12*th^3 + 8.97e-8*th^2 -5.42e-4*th + 2.02;
    end

    %Ionic Xenon Partition Function (IXPF)
    function f_p = IXPF(T)
        th = T/1; %Dimensionless Temperature T/To
        f_p = 5.8e-17*th^4 -3.71e-12*th^3 + 8.0e-8*th^2 -6.37e-4*th + 6.97;
    end

    %SAHA Equation
    function alpha = SAHA(P,T)
        %m = 2.18e-25; %Mass of Xenon [kg]
        m = 4.48e-26;     % Mass of Electron [kg]
        k = 1.38e-23; % Boltzmann's constant [J/K]
        h = 6.62607e-34;%Planck's Constant [m^2 kg/s]
        epsilon_i = 12.13; %Ionization Energy wrt. Atomic Ground State [eV]
        epsilon = 1.60218e-19*epsilon_i; %Ionization Energy [J]
        
        P = P*133.322; %torr to Pa [kg/ m s^2]
        syms alpha real
        term1 = ((2*((2*pi*m)^(3/2))*((k*T)^(5/2)))/(P*h^3));
        term2 = (IXPF(T)/AXPF(T));
        term3 = exp(-epsilon/(k*T));
        eqn1 = ((alpha^2)/(1-alpha^2)) == term1*term2*term3;
        sol = vpa(solve(eqn1, alpha));
        alpha = sol(2);
    end
end
