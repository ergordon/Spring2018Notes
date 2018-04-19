function HW6
    %Given Parameters
    MW = 39.948*0.001;             % Argon Molecular Weight [amu]
    ei = 15.76;              % Ionization Energy [eV] 
    G = 2e11;                % Scaling Constant [K^{-3/4}*m^{-3/2}]
    T_low = 5000;            % Lower Temperature Bound [K]   
    T_high = 10000;          % Higher Temperature Bound [K]
    k = 8.6173303e-5;        % Boltzmann's constant [eV/K]
    k2 = 1.38064852e-23;     % Boltzmann's Constant [J/K] or [m^2 kg s^-2 K^-1]
    %Conversions 
    
    %STP Conditions
    T_STP = 273.15; % STP Temperature [K]
    P_STP = 101325;   % STP Pressure [Pa] = [kg/ms^2]
    
    no = P_STP/(k2*T_STP); %Total Density of Heavy Particles
    %10% from STP
    no = no*.1
    
    T = T_low:1:T_high;
    alpha = zeros(length(T),1);
    eta = zeros(length(T),1);
    
    for i=1:length(T)
        %Small SAHA Equation
        alpha(i) = G*(no^(-1/2))*(T(i)^(3/4))*exp(-ei/(2*k*T(i)));
        
        %Collision Cross-Sections
        Q_en = 10e-20;                   % Electron-Neutral: Collision Cross Section [m^2]
        Q_ei = 6.5e-17/((3/2)*k*T(i))^2; % Electron-Ion: Collision Cross Section [m^2]
            
        %Collision Frequencies
        np = no*alpha(i);
        na = no*(1-alpha(i));
        
        nu_en = na*Q_en;
        nu_ei = np*Q_ei;
        
        if T(i) == 8313
            nu_ei/nu_en
            alpha(i)
            Q_ei
            np
            vv = sqrt((8*k2*T(i))/(pi*9.1e-31))  
            
            (1.6e-19*.1)/(9.1e-31*np*Q_ei*vv)
        end
        
        eta(i) = nu_ei/nu_en;
    end
    
    subplot(2,1,1)
    plot(T,alpha,'linewidth',2)
    grid on
    xlabel('Temperature [K]')
    ylabel('\alpha','FontSize',18)
    title('Degree of Ionization')
    
    subplot(2,1,2)    
    plot(T,eta,'linewidth',2)
    grid on
    xlabel('Temperature [K]')
    ylabel('\nu_{ei}/\nu_{en}','FontSize',18)
    title('Collision Frequency Ratio')
    
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','Problem2.pdf');
    
end
