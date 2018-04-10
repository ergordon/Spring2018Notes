function HW6P3
    %Given Parameters
    MW = 208.98;             % Bismuth Molecular Weight [amu]
    ei = 7.29;               % Ionization Energy [eV] 
    u = 1.2566370614e-6;     % N/A^2
    mdot_3 = 3;              % Mass Flow Rate [g/s]
    mdot_6 = 6;              % Mass Flow Rate [g/s]
    J_low = 3;               % Lower Current [kA]
    J_high = 25;             % Higher Current [kA]

    %Dimensions for the Princeton Benchmark Thruster [cm]
    PBT = [0.95, 5.1, 9.3, 6.4, 0.95, 10]; 
    %     [r_c, r_a, r_ao, r_ch, t_a, l_c]

    %Conversions 
    MW = MW*(1.66054e-27); % [amu] to [kg]
    ei = ei*(1.60218e-19); % [eV] to [Joules]
    mdot_3 = mdot_3*1e-3;  % [g/s] to [kg/s]
    mdot_6 = mdot_6*1e-3;  % [g/s] to [kg/s]
    J_low = J_low*1000;    % [kA] to [A]
    J_high = J_high*1000;  % [kA] to [A]

    J = J_low:1:J_high;

    Ct_3 = zeros(length(J),1);
    Ct_6 = zeros(length(J),1);
    T_3 = zeros(length(J),1);
    T_6 = zeros(length(J),1);
    ue_3 = zeros(length(J),1);
    ue_6 = zeros(length(J),1);
    Isp_3 = zeros(length(J),1);
    Isp_6 = zeros(length(J),1);
    jet_3 = zeros(length(J),1);
    jet_6 = zeros(length(J),1);

    ra = PBT(2);
    rc = PBT(1);

    for i=1:length(J)
        xi_3 = xi(J(i),u,ra,rc,mdot_3,ei, MW);
        xi_6 = xi(J(i),u,ra,rc,mdot_6,ei, MW);

        Ct_3(i) = C_T(mdot_3,xi_3,ra,rc);
        Ct_6(i) = C_T(mdot_6,xi_6,ra,rc);

        T_3(i) = T(Ct_3(i),u,J(i));
        T_6(i) = T(Ct_6(i),u,J(i));

        ue_3(i) = T_3(i)/mdot_3;
        ue_6(i) = T_6(i)/mdot_6;

        Isp_3(i) = ue_3(i)/9.81;
        Isp_6(i) = ue_6(i)/9.81;

        jet_3(i) = (.5*T_3(i)*ue_3(i))/1000;
        jet_6(i) = (.5*T_6(i)*ue_6(i))/1000;
    end

    JJ = linspace(3,25,length(J)); %convenient plotting for current

    figure(1);
    title('Bismuth Gas')
    clf;

    subplot(2,2,1)
    plot(JJ,Ct_3,JJ,Ct_6,'linewidth',2)
    grid on
    xlabel 'Current [kA]'
    ylabel 'Thrust Coefficient (C_T)'
    title 'Thrust Coefficient'
    xlim([3,25])
    ylim([0,8])
    leg1 = legend('3 g/s','6 g/s');

    subplot(2,2,2)
    plot(JJ,T_3,JJ,T_6,'linewidth',2)
    grid on
    xlabel 'Current [kA]'
    ylabel 'Thrust [N]'
    title 'Thrust'
    xlim([3,25])
    leg2 = legend('3 g/s','6 g/s','Location','northwest')

    subplot(2,2,3)
    plot(JJ,Isp_3,JJ,Isp_6,'linewidth',2)
    grid on
    xlabel 'Current [kA]'
    ylabel 'I_{SP} [sec]'
    title 'I_{SP}'
    xlim([3,25])
    leg3 = legend('3 g/s','6 g/s','Location','northwest')

    subplot(2,2,4)
    plot(JJ,jet_3,JJ,jet_6,'linewidth',2)
    grid on
    xlabel 'Current [kA]'
    ylabel 'Jet Power [kW]'
    title 'Jet Power'
    xlim([3,25])
    leg4 = legend('3 g/s','6 g/s','Location','northwest')
    
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','Problem3.pdf');
end

function nu = nu(mdot)
    mdot_star = 0.066;
    nu = mdot/mdot_star;
end

function xi = xi(J,u,ra,rc,mdot,ei, MW)
    xi = (J*sqrt((u/(4*pi))*log(ra/rc)))/(sqrt(mdot)*(((2*ei)/MW)^(1/4)));
end

function C_T = C_T(mdot,xi,ra,rc)
    C_T = (nu(mdot)/(xi^4)) + log((ra/rc) + (xi^2));
end

function T = T(C_T,u,J)
    T = C_T*(u/(4*pi))*(J^2);
end