clc; clear;

% Initialize variables
m = 4.48e-26;     % Mass of Electron [kg]
k = 1.38e-23;     % Boltzmann's constant [J/K]
Te = [5; 50];           % Particle Temperature [eV]




% Convert eV to K
T = 1.16045221e4*Te; %Temperature [K]

% Compute two frequently used constants
C1 = [4*pi*(m/(2*pi*k*T(1)))^(3/2);4*pi*(m/(2*pi*k*T(2)))^(3/2);]
C2 = [m/(2*k*T(1));m/(2*k*T(2))];
C3 = [(m/(2*pi*k*T(1)))^(1/2);(m/(2*pi*k*T(2)))^(1/2)];

c_mp1 = sqrt((2*k*T(1))/(m)) %Most Probable Speed
c_mp2 = sqrt((2*k*T(2))/(m))

c = [-3*c_mp2:50:3*c_mp2];  % Particle Velocity Magnitude [m/s]
for i=1:length(c)
    %Maxwellian Speed Distribution Function
    if c(i)<0
        chiM(1,i)=0;
        chiM(2,i)=0;
    else
        chiM(1,i) = C1(1) * c(i)^2 * exp(-C2(1) * c(i)^2);
        chiM(2,i) = C1(2) * c(i)^2 * exp(-C2(2) * c(i)^2);
    end
    
    %Maxwellian Velocity Distribution Function
    fM(1,i) = C3(1) * exp(-C2(1) * c(i)^2);
    fM(2,i) = C3(2) * exp(-C2(2) * c(i)^2);
    %fM1 = C3 * exp(-C2 * c(i)^2);
    %fM(i) = norm([fM1, fM1, fM1]);
end
hold on
plot(c/c_mp1,chiM(1,:)*c_mp1,'-b','LineWidth',2) %normalize with cmp.
plot(c/c_mp1,chiM(2,:)*c_mp1,'-r','LineWidth',2) %normalize with cmp.
plot(c/c_mp1,fM(1,:)*c_mp1,'--b','LineWidth',2)
plot(c/c_mp1,fM(2,:)*c_mp1,'--r','LineWidth',2)
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
    'Maxwellian Speed Distribution T_i = 50 eV',...
    'Maxwellian Velocity Distribution T_i = 5 eV',...
    'Maxwellian Velocity Distribution T_i = 50 eV',...
    'Location','northwest')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','electronBoth.pdf');
