function hw6p3

clc;

J1 = 12;
J2 = 14;
J3 = 8; 
Jw = 1; 
n = .0011;


A = [      0                 0               1               0          0     0;...
           0                 0               0               1          0     0;...
     -4*n^2*(J2-J3)/J1       0               0         n*(1-(J2-J3)/J1) 0     0;...
           0           n^2*(J1-J2)/J3 -n*(1+(J1-J2)/J3)      0          0     0;...
           0                 0               0               0          0     n;...
           0                 0               0               0         -n     0];
       
B = [0           0;...
     0           0;...
     1/J1        0;...
     0        1/J3;...
     -1/Jw       0;...
     0       -1/Jw];
 
rank(ctrb(A,B)); %6 is full rank, this full rank 
x0 = [0.1 0.5 0 0 0 0]';

Q = eye(6);
R = 10^0.*eye(2); 

t = linspace(0,100,100);

% Part B
[k,s,e] = lqr(A,B,Q,R)


for i=1:length(t)
    x(:,i) = expm((A-B*k)*t(i))*x0;
end

figure(1) %x(t)
plot(t,x,'linewidth',2);
hold on;
legend('\theta_1','\theta_3','\theta_1dot','\theta_3dot','v_1','v_3')
title('x(t) V. Time');
set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','HW6P3C1.pdf')
    
figure(2)
subplot(3,2,1)
plot(t,x(1,:),'linewidth',2);
title('\theta_1 V. Time');
legend('x_1(t) - \theta_1');

subplot(3,2,3)
plot(t,x(2,:),'linewidth',2);
title('\theta_3 V. Time');
legend('x_2(t) - \theta_3');

subplot(3,2,2)
plot(t,x(3,:),'linewidth',2);
title('\theta_1 dot V. Time');
legend('x_3(t) - \theta_1 dot ');

subplot(3,2,4)
plot(t,x(4,:),'linewidth',2);
title('\theta_3 dot V. Time');
legend('x_3(t) - \theta_3 dot');

subplot(3,2,5)
plot(t,x(5,:),'linewidth',2);
title('v_1 V. Time');
legend('x_5(t) - v_1');

subplot(3,2,6)
plot(t,x(6,:),'linewidth',2);
title('v_3 V. Time');
legend('x_6(t) - v_3');

set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','HW6P3C2.pdf')