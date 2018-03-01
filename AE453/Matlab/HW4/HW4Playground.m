function HW4Playground

    % Initialize variables
    m = 4.48e-26;     % Mass of Electron [kg]
    k = 1.38e-23;     % Boltzmann's constant [J/K]
    temp = 50;
    Te = [5; temp];           % Particle Temperature [eV]

    % Convert eV to K
    T = 1.16045221e4*Te; %Temperature [K]
    
    f = figure('Visible','off');
    
    [x,y1,y2]=pf(temp)
    chiM = y1;
    fM = y2;
    hold on
    plot(c/c_mp1,chiM(1,:)*c_mp1,'-b','LineWidth',2) %normalize with cmp.
    plot(c/c_mp1,chiM(2,:)*c_mp1,'-r','LineWidth',2) %normalize with cmp.
    plot(c/c_mp1,fM(1,:)*c_mp1,'--b','LineWidth',2)
    plot(c/c_mp1,fM(2,:)*c_mp1,'--r','LineWidth',2)
    xL = xlim;
    yL = ylim;
    line(xL, [0 0],'color','k','linewidth',2) %x-axis 
    line([0 0], yL,'color','k','linewidth',2) %y-axis
    hold off
    title('Maxwellian Distribution Functions for Electron', 'FontSize', 18);
    ylabel('Probability', 'FontSize', 13)
    xlabel('c_1/(2kT/m)^{1/2} or C/(2kT/m)^{1/2}', 'FontSize', 13)
    grid on
    grid minor
    legend('Maxwellian Speed Distribution T_i = 5 eV',...
        'Maxwellian Speed Distribution T_i = 50 eV',...
        'Maxwellian Velocity Distribution T_i = 5 eV',...
        'Maxwellian Velocity Distribution T_i = 50 eV',...
        'Location','northwest')

 % Create slider
    sld1 = uicontrol('Style', 'slider',...
        'Min',5,'Max',100,'Value',50,...
        'Position', [200 0 200 20],...
        'Callback', @Templim); 

    % Add a text uicontrol to label the slider.
    txt1 = uicontrol('Style','text',...
        'Position',[100 2 100 20],...
        'String','Temperature [5:100]');

    % Make figure visble after adding all components
    f.Visible = 'on';

    function plotgraph(source,event)
        temp = source.Value;
    end

    function Templim(source,event)
        temp = source.Value;
        [x,y1,y2]=pf(temp)
        chiM = y1;
        fM = y2;
        hold off
        plot(c/c_mp1,chiM(1,:)*c_mp1,'-b','LineWidth',2)
        hold on
        plot(c/c_mp1,chiM(2,:)*c_mp1,'-r','LineWidth',2) %normalize with cmp.
        plot(c/c_mp1,fM(1,:)*c_mp1,'--b','LineWidth',2)
        plot(c/c_mp1,fM(2,:)*c_mp1,'--r','LineWidth',2)
        hold off
        xL = xlim;
        yL = ylim;
        line(xL, [0 0],'color','k','linewidth',2) %x-axis 
        line([0 0], yL,'color','k','linewidth',2) %y-axis
        title('Maxwellian Distribution Functions for Electron', 'FontSize', 18);
        ylabel('Probability', 'FontSize', 13)
        xlabel('c_1/(2kT/m)^{1/2} or C/(2kT/m)^{1/2}', 'FontSize', 13)
        grid on
        grid minor
        legend('Maxwellian Speed Distribution T_i = 5 eV',...
            sprintf('Maxwellian Speed Distribution T_i = %d eV',round(temp)),...
            'Maxwellian Velocity Distribution T_i = 5 eV',...
            sprintf('Maxwellian Velocity Distribution T_i = %d eV',round(temp)),...
            'Location','northwest')
    
        txt1 = uicontrol('Style','text',...
        'Position',[100 2 100 20],...
        'String',['Temperature=' num2str(temp)]);
    end

    function [x,y1,y2] = pf(temp)
        % Initialize variables
        m = 4.48e-26;     % Mass of Electron [kg]
        k = 1.38e-23;     % Boltzmann's constant [J/K]
        Te = [5; temp];           % Particle Temperature [eV]
        T = 1.16045221e4*Te; %Temperature [K]
        C1 = [4*pi*(m/(2*pi*k*T(1)))^(3/2);4*pi*(m/(2*pi*k*T(2)))^(3/2)];
        C2 = [m/(2*k*T(1));m/(2*k*T(2))];
        C3 = [(m/(2*pi*k*T(1)))^(1/2);(m/(2*pi*k*T(2)))^(1/2)];

        c_mp1 = sqrt((2*k*T(1))/(m)); %Most Probable Speed
        c_mp2 = 2.6739e+04;

        c = [-3*c_mp2:50:3*c_mp2];  % Particle Velocity Magnitude [m/s]
        for i=1:length(c)
            %Maxwellian Speed Distribution Function
            if c(i)<0
                chiMM(1,i)=0;
                chiMM(2,i)=0;
            else
                chiMM(1,i) = C1(1) * c(i)^2 * exp(-C2(1) * c(i)^2);
                chiMM(2,i) = C1(2) * c(i)^2 * exp(-C2(2) * c(i)^2);
            end

            %Maxwellian Velocity Distribution Function
            fMM(1,i) = C3(1) * exp(-C2(1) * c(i)^2);
            fMM(2,i) = C3(2) * exp(-C2(2) * c(i)^2);
        end
        x = c;
        y1 = chiMM;
        y2 = fMM;
    end
end
