function hw4p1
t = [0:1:250]/1000;
i=1;
dzdt = 0;
z = 0;
V = 0;

while i<1+length(t)
    if i<100
        dzdt = [dzdt 2];
        z = [z zloc(dzdt(i),t(i))];
        V = [V Vz(dzdt(i),z(i))];
    elseif (i>100)&(i<151)
        dzdt = [dzdt 0];
        z = [z zloc(dzdt(100),t(100))+zloc(dzdt(i),t(i))];
        V = [V Vz(dzdt(i),z(i))];
    elseif i>150
        dzdt = [dzdt -2];
        z(150)
        zloc(dzdt(i),t(i))
        z = [z z(150)+zloc(dzdt(i),(t(i)-.15))];
        V = [V Vz(dzdt(i),z(i))];
    end
    i=i+1;
end
figure(3)
plot(t,z)
figure(1)
plot(t*1000,V,'LineWidth',2)
title('Voltage Trace V. Time with Varying Distance', 'FontSize', 18);
ylabel('Voltage [ V ]', 'FontSize', 13)
xlabel('Time [ms]', 'FontSize', 13)
grid on
grid minor

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem1VoltageVTime.pdf');

AxialDistance = [0:.0002:.200]
for j=1:length(AxialDistance)
        VV(1,j)=Vz(2,AxialDistance(j));
        VV(2,j)=Vz(-2,AxialDistance(j));
end
figure(2)
plot(AxialDistance*1000,VV(1,:),'-b','LineWidth',2)
hold on
plot(AxialDistance*1000,VV(2,:),'-r','LineWidth',2)
title('Voltage Trace V. Axial Distance', 'FontSize', 18);
ylabel('Voltage [ V ]', 'FontSize', 13)
xlabel('Axial Distance in z [mm]', 'FontSize', 13)
grid on
grid minor
legend('Table Velocity = 2[m/s]',...
    'Table Velocity = -2[m/s]',...
    'Location','northeast')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','Problem1VoltageVDistance.pdf');

end
function V = Vz(dzdt,z)
    n = 50;
    A = 1.963e-5;
    mu = 1.256e-6;
    JN = 100;
    a = 0.05;
    Constant1 = ((n*A*3*mu*JN*a^2)/2)*dzdt;
    V = Constant1*(z/((a^2+z^2)^(5/2)));
end
function z = zloc(dzdt,t)
    z = dzdt*t;
end