function hw4
clf;
%hw4a(0.044,500) %Part A
hw4b(0.044,500) %Part B
%hw4a(10,500) %Part C
hw4b(10,500) %Part C
hw4b(0.044,10000) %Part B

function hw4a(J3,tMax)
format long

R0 = [1 0 0;0 1 0; 0 0 1]; % Axis Initially Aligned
woz = rpm2rads(0.01);      % [rad/s]
w0 = [0;0;woz];            % Initial Spin Conditions [rad/s]
J = [4.83 0 0; 0 4.83 0; 0 0 J3]; % Principle Moment of Inertia [kg m^2]
Jd = 1;            % Moment of Inertia [kg m^2]
n = rpm2rads(100); % [rad/s]
c = 1;             % Coefficient of Friction on Wheel

% Define the Time Interval
tRate = 1;
t = linspace(0,tMax,tMax*tRate);
t1 = linspace(0,tMax,8886);

% Integrate Angular Rates and Euler's Equations
[t,x] = ode45(@(t,x) fa(t,x,J,Jd,n,c,t1),t,[RtoX(R0); w0]);  

figure
start = 1
plot(t(start:end),x(start:end,10),t(start:end),x(start:end,11),t(start:end),x(start:end,12)); 
title('Angular Velocity V. Time')
ylabel('w (rad/s)')
xlabel('time (s)')
legend('w1','w2','v')
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf',sprintf('LinearEquationsJ3%dtMax%d.pdf',J3,tMax));


function hw4b(J3,tMax)
woz = rpm2rads(.01); % [rad/s]
n = rpm2rads(100);   % [rad/s]
x0 = [0;0;woz;n];    % Initial Spin Conditions [w1 w2 w3 v]
J = [4.83 0 0; 0 4.83 0; 0 0 J3]; % Principle Moment of Inertia of Craft [kg m^2]
Jd = 1;    % Momet of Inertia of Wheel [kg m^2]
c = 1;     % Coefficient of Friction [N m /rad/s]


% Define the Time Interval
tRate = 1;
t = linspace(0,tMax,tMax*tRate);

% Integrate Angular Rates and Euler's Equations
[t,x] = ode45(@(t,x) fb(t,x,J,Jd,c),t,x0);  

figure 
subplot(2,2,1)
plot(t,x(:,1)); 
title('Angular Velocity (w_1) V. Time')
ylabel('w_1 (rad/s)')
xlabel('time (s)')

subplot(2,2,2)
plot(t,x(:,2)); 
title('Angular Velocity (w_2) V. Time')
ylabel('w_2 (rad/s)')
xlabel('time (s)')

subplot(2,2,3)
plot(t,x(:,3)); 
title('Angular Velocity (w_3) V. Time')
ylabel('w_3 (rad/s)')
xlabel('time (s)')

subplot(2,2,4)
plot(t,x(:,4)); 
title('Angular Velocity (v) V. Time')
ylabel('v (rad/s)')
xlabel('time (s)')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf',sprintf('NonLinearEquationsJ3%dtMax%d.pdf',J3,tMax));


function xdot = fa(t,x,J,Jd,n,c,t1)
R = XtoR(x(1:9,1));
w = x(10:end,1);

Jt = J(1,1);
Ja = J(3,3);
Jd = Jd;
n = n; 

Rdot = R*skew(w,n);
A = [0, -n*(Ja-Jt)/Jt, n*Jd/Jt;...
    n*(Ja-Jt)/(Jt-Jd), 0, c/(Jt-Jd);...
    -n*(Ja-Jt)/(Jt-Jd), 0, -c*Jt/(Jd*(Jt-Jd))];
w = [w(1); w(2); w(3)];
wdot = A*w;

xdot = [RtoX(Rdot); wdot];

function xdot = fb(t,x,J,Jd,c)
w = x; % [w1 w2 w3 v]
Jt = J(1,1);
Ja = J(3,3);
n = w(3);
v = w(4);

wdot = [((Ja-Jt)*w(2)*w(3)-Jd*w(3)*v)/-Jt; ...
    (1/(Jt-Jd))*(w(1)*w(3)*(Ja-Jt)+c*v);...
    (-Jd*w(1)*v)/Ja;...
    (-c*v/Jd)-((1/(Jt-Jd))*(w(1)*w(3)*(Ja-Jt)+c*v))];

 xdot = wdot;

function S = skew(w,n)
S = [0 -n w(2); n 0 -w(1); -w(2) w(1) 0];

function X = RtoX(R)
X = reshape(R,9,1);

function R = XtoR(X)
R = reshape(X,3,3);

function rad = rpm2rads(rpm)
rad = rpm*0.104719755;