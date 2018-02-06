
clear
clc

alpha=200
n=0.80
tt=6.307e+7;
vc2 = 2*tt*n*alpha/(1000^2);
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
plot(c,mplmo1,'r-',c,mplmo2,'b-',c,mplmo3,'g-',c,mplmo4,'k-',c,mplmo5,'m-');
legend('dV1 = 1 km/s','dV2 = 5 km/s','dV3 = 10 km/s','dV4 = 20 km/s','dV5 = 50 km/s');
xlim([10,1000]);
ylim([0,1]);
xlabel 'c: m/s'
ylabel ' mpl/mo'
title 'Payload Ratio Vs Effective Exhaust Velocity'