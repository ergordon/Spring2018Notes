clc; clear;

% Initialize variables
m = 4.48e-26;     % Mass of nitrogen molecule [kg]
k = 1.38e-23;     % Boltzmann's constant [J/K]
Te = .5;           % Particle Temperature [eV]



% Convert eV to K
T = 1.16045221e4*Te; %Temperature [K]

% Compute two frequently used constants
C1 = 4*pi*(m/(2*pi*k*T))^(3/2);
C2 = m/(2*k*T);
C3 = (m/(2*pi*k*T))^(1/2);

c_mp = sqrt((2*k*T)/(m)); %Most Probable Speed
c = [-3*c_mp:100:3*c_mp];  % Particle Velocity Magnitude [m/s]
for i=1:length(c)
    %Maxwellian Speed Distribution Function
    if c(i)<0
        chiM(i)=0;
    else
        chiM(i) = C1 * c(i)^2 * exp(-C2 * c(i)^2);
    end
    
    %Maxwellian Velocity Distribution Function
    fM(i) = C3 * exp(-C2 * c(i)^2);
    %fM1 = C3 * exp(-C2 * c(i)^2);
    %fM(i) = norm([fM1, fM1, fM1]);
end
hold on
plot(c/c_mp,chiM*c_mp,'LineWidth',2) %normalize with cmp.
plot(c/c_mp,fM*c_mp,'LineWidth',2)

xL = xlim;
yL = ylim;
line(xL, [0 0],'color','k','linewidth',2) %x-axis 
line([0 0], yL,'color','k','linewidth',2) %y-axis

hold off

title(sprintf('Maxwellian Distribution Functions for Xenon Ions with Ti=%d [eV]',Te), 'FontSize', 18)
ylabel('Probability', 'FontSize', 13)
xlabel('c_1/(2kT/m)^{1/2} or C/(2kT/m)^{1/2}', 'FontSize', 13)
grid on
grid minor
legend('Maxwellian Speed Distribution \chi(c)','Maxwellian Velocity Distribution \phi(c_1)','Location','northwest')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf',sprintf('distXe%deV.pdf',Te));
