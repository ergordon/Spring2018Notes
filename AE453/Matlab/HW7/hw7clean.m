function hw7clean
clc; clear;
l = 0.05; % Channel Length [m]

%Constant Density
tMaxCD = 4.69e-6;
tspan = linspace(0, tMaxCD,100);
[tCD,xx] = ode23s(@snowplowi,tspan,[0.00001,0]);
xCD = xx(:,1);
xdotCD = xx(:,2);

%Linearly Decreasing Density
tMaxLDD = 3.8311e-6;
tspan = linspace(0, tMaxLDD, 100);
[tLDD,xx] = ode23s(@snowplowii,tspan,[0.00001,0]);
xLDD = xx(:,1);
xdotLDD = xx(:,2);

%Square Root Dependence,Decreasing Density 
tMaxSQD = 4.19633e-6;
tspan = linspace(0, tMaxSQD,100);
[tSQD,xx] = ode23s(@snowplowiii,tspan,[0.00001;0]);

xSQD = xx(:,1);
xdotSQD = xx(:,2);

figure(1)
subplot(1,2,1)
plot(tCD,xCD,'linewidth',2)
title(sprintf('Position V. Time \n Constant Density Snowplow Model'))
xlabel 'Time     [sec]'
ylabel 'Position [m]'
grid on
grid minor
subplot(1,2,2)
plot(tCD,xdotCD,'linewidth',2)
title(sprintf('Velocity V. Time \n Constant Density Snowplow Model'))
xlabel('Time     [sec]')
ylabel 'Velocity [m/s]'
grid on
grid minor
hold off
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartCi.pdf')

figure(2)
subplot(1,2,1)
plot(tLDD,xLDD,'linewidth',2)
title(sprintf('Position V. Time \n Linear Decreasing Density Snowplow Model'))
xlabel 'Time     [sec]'
ylabel 'Position [m]'
grid on
grid minor
subplot(1,2,2)
plot(tLDD,xdotLDD,'linewidth',2)
title(sprintf('Velocity V. Time \n  Linear Decreasing Density Snowplow Model'))
xlabel 'Time     [sec]'
ylabel 'Velocity [m/s]'
grid on
grid minor
hold off
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartCii.pdf')

figure(3)
subplot(1,2,1)
plot(tSQD,xSQD,'linewidth',2)
title(sprintf('Position V. Time \n Square Root Dependence,Decreasing Density Snowplow Model'))
xlabel 'Time     [sec]'
ylabel 'Position [m]'
grid on
grid minor
subplot(1,2,2)
plot(tSQD,xdotSQD,'linewidth',2)
title(sprintf('Velocity V. Time \n  Square Root Dependence,Decreasing Density Snowplow Model'))
xlabel 'Time     [sec]'
ylabel 'Velocity [m/s]'
grid on
grid minor
hold off
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartCiii.pdf')

figure(4)
plot(tCD,xCD,'linewidth',2)
hold on
plot(tLDD,xLDD,'linewidth',2)
hold on
plot(tSQD,xSQD,'linewidth',2)
title(sprintf('Position V. Time \n Snowplow Model'))
legend('Constant Density','Linear Decreasing Density','Square Root Dependence,Decreasing Density','location','bestoutside')
xlabel 'Time     [sec]'
ylabel 'Position [m]'
grid on
grid minor
hold off
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartD.pdf')

figure(5)
plot(tCD,xdotCD,'linewidth',2)
hold on
plot(tLDD,xdotLDD,'linewidth',2)
hold on
plot(tSQD,xdotSQD,'linewidth',2)
title(sprintf('Velocity V. Time \n Snowplow Model'))
legend('Constant Density','Linear Decreasing Density','Square Root Dependence,Decreasing Density','location','bestoutside')
xlabel 'Time     [sec]'
ylabel 'Velocity [m/s]'
grid on
grid minor
hold off
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem2PartE.pdf')

end


function xdot = snowplowi(t,x)
l = 0.05;                   % Channel Length [m]
area = 0.03*0.03;           % Area of Channel
R_A = 208.13;               % Specific Gas Constant for Argon [J/kgK]
P_o = 66.661185;            % Peak Pressure [Pa]
T_o = 273;                  % Peak Temperature [K]
rho = P_o / (R_A * T_o);    % Total Density of Heavy Particles
J = 20e3;                   % Initial Current Sheet Disharge [A]
Lprime = 0.6E-6;            % Distributed Inductance [H]
F = 0.5*Lprime*J^2;         % Thrusting Force
m = area*rho*x(1);
mdot = area*rho*x(2);
xdot = [x(2); (F-(mdot*x(2)))/m];
end


function xdot = snowplowii(t,x)
l = 0.05;                   % Channel Length [m]
area = 0.03*0.03;           % Area of Channel
R_A = 208.13;               % Specific Gas Constant for Argon [J/kgK]
P_o = 66.661185;            % Peak Pressure [Pa]
T_o = 273;                  % Peak Temperature [K]
rho = P_o / (R_A * T_o);    % Total Density of Heavy Particles
J = 20e3;                   % Initial Current Sheet Disharge [A]
Lprime = 0.6E-6;            % Distributed Inductance [H]
F = 0.5*Lprime*J^2;         % Thrusting Force
m = area*rho*(x(1)-(((x(1))^2)/(2*l)));
mdot = area*rho*x(2)*(1-(x(1)/l));
xddot = F\[mdot m]
xdot = [x(2); (F-(mdot*x(2)))/m]
end

function xdot = snowplowiii(t,x)
l = 0.05;                   % Channel Length [m]
area = 0.03*0.03;           % Area of Channel
R_A = 208.13;               % Specific Gas Constant for Argon [J/kgK]
P_o = 66.661185;            % Peak Pressure [Pa]
T_o = 273;                  % Peak Temperature [K]
rho = P_o / (R_A * T_o);    % Total Density of Heavy Particles
J = 20e3;                   % Initial Current Sheet Disharge [A]
Lprime = 0.6E-6;            % Distributed Inductance [H]
F = 0.5*Lprime*J^2;         % Thrusting Force
xx = x(1);
m = area*rho*(-(2/3)*(l - x(1))*sqrt(1 - x(1)/l)+((2/3)*l));
mdot = area*rho*sqrt(1-(x(1)/l))*x(2);
xdot = [x(2); (F-(mdot*x(2)))/m];
end