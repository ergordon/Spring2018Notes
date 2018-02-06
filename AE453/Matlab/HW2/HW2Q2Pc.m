clc; clear; 

alpha=200 % Specific Power [kW/kg] 
n=0.80 % Thruster Efficiency [%]
tt=6.307e+7; % Thrusting Time [sec]
vc2 = 2*tt*n*alpha/(1000^2); % Charachteristic [km/s]
c = linspace(10,1000,1000); 

dv1 = 1;
dv2 = 5;
dv3 = 10;
dv4 = 20
dv5 = 50;

mplmo1 = zeros(length(c),1);
mplmo2 = zeros(length(c),1);
mplmo3 = zeros(length(c),1);
mplmo4 = zeros(length(c),1);
mplmo5 = zeros(length(c),1);

for i = 1:length(c)
    mplmo1(i) = exp(-dv1/c(i))-(((c(i)^2)/(vc2))*(1-exp(-dv1/c(i))));    
    mplmo2(i) = exp(-dv2/c(i))-(((c(i)^2)/(vc2))*(1-exp(-dv2/c(i))));
    mplmo3(i) = exp(-dv3/c(i))-(((c(i)^2)/(vc2))*(1-exp(-dv3/c(i))));
    mplmo4(i) = exp(-dv4/c(i))-(((c(i)^2)/(vc2))*(1-exp(-dv4/c(i))));
    mplmo5(i) = exp(-dv5/c(i))-(((c(i)^2)/(vc2))*(1-exp(-dv5/c(i))));
end

p = plot(c,mplmo1,...
         c,mplmo2,...
         c,mplmo3,...
         c,mplmo4,...
         c,mplmo5);
     p(1).LineWidth = 2;
     p(2).LineWidth = 2;
     p(3).LineWidth = 2;
     p(4).LineWidth = 2;
     p(5).LineWidth = 2;
     title('Payload Ratio v. Effective Exhaust Velocity', 'FontSize', 24)
     ylabel('Payload Ratio [kg/kg]', 'FontSize', 18)
     xlabel('Effective Exhaust Velocity [km/s]', 'FontSize', 18)
     xlim([10,1000]);
     ylim([0,1]);
     legend('Delta-V = 1 [km/s]','Delta-V = 5 [km/s]','Delta-V = 10 [km/s]', 'Delta-V = 20 [km/s]', 'Delta-V = 50 [km/s]')
     set(gcf,'paperorientation','landscape');
     set(gcf,'paperunits','normalized');
     set(gcf,'paperposition',[0 0 1 1]);
     print(gcf,'-dpdf','graph.pdf');

%}