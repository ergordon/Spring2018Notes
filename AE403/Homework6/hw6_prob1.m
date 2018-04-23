function hw6_prob1
clear; clc;
%% Problem 1 Part D
wn = [10^-2 10^-1 10^0 10^1 10^2];
J = 10;
Jw = 1; 
A = [0 1; 0 0];
B = [0 1/J]';
x0 = [0.1 0]';
v0 = 0;

for j=1:length(wn)

    tRate = 30;
    tMax = 1;
    t = linspace(0,tMax,tMax*tRate);

    k = [wn(j)^2*J 2*J/sqrt(2)*wn(j)];
    
    for i=1:length(t)
        x(:,i) = expm((A-B*k)*t(i))*x0;
    end
    
    for i=1:length(t)
        v(:,i) = -k*(1/Jw)*inv(A-B*k)*(expm((A-B*k)*t(i))-eye(2))*x0;
    end
    
    theta(:,j) = x(1,:);
    uu(:,j)=-k*x;
    vv(:,j)=v;
end
    figure(1)
    subplot(3,1,1)
    plot(t,theta,'linewidth',2)
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('\theta(t) - Problem 1 Part D')
    ylabel('\theta')
    xlabel('Time')
    
    subplot(3,1,2)
    plot(t,uu,'linewidth',2)
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('u(t) - Problem 1 Part D')
    ylabel('u')
    xlabel('Time')
    
    subplot(3,1,3)
    plot(t,vv,'linewidth',2)
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('v(t) - Problem 1 Part D')
    ylabel('v')
    xlabel('Time')
    
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','HW6P1D.pdf')

%% Problem 1 Part E
R = [10^-2 10^-1 10^0 10^1 10^2];
J = 10;
Jw = 1;
A = [0 1; 0 0];
B = [0 1/J]'; 
x0 = [0.1 0]';
v0 = 0;
Q = eye(2);

for j=1:length(R)

    t = linspace(0,10,1400);
    k = lqr(A,B,Q,R(j));

    for i=1:length(t)
        x(:,i) = expm((A-B*k)*t(i))*x0;
    end
    
    for i=1:length(t)
        v(:,i) = -k*(1/Jw)*inv(A-B*k)*(expm((A-B*k)*t(i))-eye(2))*x0;
    end
    
    thetae(:,j) = x(1,:);
    uue(:,j)=-k*x;
    vve(:,j)=v;
end 
    figure(2)
    subplot(3,1,1)
    plot(t,thetae,'linewidth',2)
    hold on;
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('\theta(t) - Problem 1 Part E')
    ylabel('\theta')
    xlabel('Time')

    subplot(3,1,2)
    plot(t,uue,'linewidth',2)
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('u(t) - Problem 1 Part E')
    ylabel('u')
    xlabel('Time')

    subplot(3,1,3)
    plot(t,vve)
    hold on;
    legend('\omega_n = 10e-2','\omega_n = 10e-1','\omega_n = 10e0','\omega_n = 10e1','\omega_n = 10e2','Location','bestoutside')
    title('v(t) - Problem 1 Part E')
    ylabel('v')
    xlabel('Time')
    
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','HW6P1E.pdf')
