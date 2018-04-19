function HW6P1
J = 10;
J_w = 1;

%[Theta dTheta]
x_0=[0.1; 0];

A = [0 1; 0 0];
B = [0; 1/J];


tRate = 30;
tMax = 50;
t = linspace(0,tMax,tMax*tRate);

wn = [10e-2; 10e-1; 10e0; 10e1; 10e2];
for i=1:length(wn)
    [t,x] = changeTheta(A,B,J,J_w,t,x_0,wn(i));
    theta(:,i) = x(:,1);
    K = [J*(wn(i))^2 ((2*J*wn(i))/(sqrt(2)))];
    C = A-B*K;
    v(:,i) = K*inv(C)*(expm(C*t)-eye(size(C)))*x_0;
    for j=1:length(x)
        u(j,i)=-[J*(wn(i))^2 ((2*J*wn(i))/(sqrt(2)))]*x(j,:)';
    
    end
end

figure(1)
plot(t,theta,'linewidth',2)
title('Theta V. Time', 'FontSize', 18)

legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')

figure(2)
plot(t,u,'linewidth',2)

figure(3)
plot(t,v,'linewidth',2)
%{
[t,x] = changer(A,B,J,J_w,t,x_0,1);
theta = x(:,1);
plot(t,theta)
%}

function [t, state] = changeTheta(A,B,J,J_w,t,x_0,wn)
    K = [J*wn.^2 ((2*J*wn)/(sqrt(2)))];
[t,state] = ode45(@(t,x) f(A,B,K,t,x_0),t,x_0);


function xdot = f(A,B,K,t,x_0)
xdot = expm((A-B*K)*t)*x_0;