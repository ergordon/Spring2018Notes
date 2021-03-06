function hw3_prob5

clc;

% %%
%
% YOUR CODE HERE TO DEFINE THE INITIAL CONDITIONS AND THE MOMENT OF INERTIA
% MATRIX AND THE TIME REQUIRED FOR SPIN-UP
%
w0x = -4
J1 = 4600; % kgm^3
J2 = 4400; % kgm^3
J3 = 750;  % kgm^3
R0 = eye(3);
w0 = [1*10^w0x; 0; 0];
J = [J1 0 0; 0 J2 0; 0 0 J3];
tMax = 7.5*2*pi;
%
% %%

% DEFINE THE TIME INTERVAL
tRate = 30;
t = linspace(0,tMax,tMax*tRate);

% INTEGRATE THE ANGULAR RATE EQUATIONS **AND** EULER'S EQUATIONS
[t,x] = ode45(@(t,x) f(t,x,J),t,[RtoX(R0); w0]);

% CREATE A BOX AND SOME AXES
Hb = 1;
Mb = 6*(J(1,1)+J(2,2)-J(3,3))/(Hb^2);
Wb = sqrt((6/Mb)*(J(1,1)-J(2,2)+J(3,3)));
Lb = sqrt((6/Mb)*(-J(1,1)+J(2,2)+J(3,3)));
[p1,faces,colors] = makebox(Lb,Wb,Hb);
pAxis = [0 0; 0 0; 0 1];
pAngVel = [zeros(3,1) w0/norm(w0)];
pAngMom = [zeros(3,1) J*w0/norm(J*w0)];

finalTime = t(end)
finalRot = reshape(x(end,1:9),3,3)
finalW = reshape(x(end,10:12),3,1)

% SETUP THE PLOT
figure(1);
clf;
axis(1.25*[-1 1 -1 1 -1 1]);
axis equal;
hold on;
plot3(0,0,0,'k.','markersize',16);
hAxis = line(pAxis(1,:),pAxis(2,:),pAxis(3,:));
set(hAxis,'linewidth',2,'color','r');
hAngVel = line(pAngVel(1,:),pAngVel(2,:),pAngVel(3,:));
set(hAngVel,'linewidth',2,'color','g');
hAngMom = line(pAngMom(1,:),pAngMom(2,:),pAngMom(3,:));
set(hAngMom,'linewidth',2,'color','b');
hBox = patch('Vertices',p1','Faces',faces,...
          'CData',colors,'FaceColor','flat');
hTitle = title(sprintf('t = %4.2f',0));
lighting flat
light('Position',[0 -2 -1])
light('Position',[0 -2 1])
xlabel('x');
ylabel('y');
zlabel('z');

% ANIMATE THE RESULTS
w = x(:,10:end);
x = x(:,1:9);
firsttime = 1;
i = 1;
dt = max(t)/(length(t)-1);
tic;
AngularMomentum = [1;0;0];
while (i<length(t))
    if (toc > dt)
        tic;
        i = i+1;
        R = XtoR(x(i,:));
        p0 = R*p1;
        set(hBox,'Vertices',p0');
        pAxis0 = R*pAxis;
        pAngVel = [zeros(3,1) w(i,:)'/norm(w(i,:)')];
        pAngVel0 = R*pAngVel;
        pAngMom = [zeros(3,1) J*(w(i,:)')/norm(J*(w(i,:)'))];
        pAngMom0 = R*pAngMom;
        set(hAxis,'xdata',pAxis0(1,:),'ydata',pAxis0(2,:),'zdata',pAxis0(3,:));
        set(hAngVel,'xdata',pAngVel0(1,:),'ydata',pAngVel0(2,:),'zdata',pAngVel0(3,:));
        set(hAngMom,'xdata',pAngMom0(1,:),'ydata',pAngMom0(2,:),'zdata',pAngMom0(3,:));
        set(hTitle,'string',sprintf('t = %4.2f',t(i)));
        ppmom = [pAngMom0(1,2); pAngMom0(2,2) ;pAngMom0(3,2)];
        norm(ppmom)
        AngularMomentum = [AngularMomentum ppmom];
        if (firsttime)
            firsttime=0;
            pause(0.5);
        end
        drawnow;
    end
end
RR1 = x(:,1);
RR2 = x(:,2);
RR3 = x(:,3);
RR4 = x(:,4);
RR5 = x(:,5);
RR6 = x(:,6);
RR7 = x(:,7);
RR8 = x(:,8);
RR9 = x(:,9);
ww1 = w(:,1);
ww2 = w(:,2);
ww3 = w(:,3);
AM1 = AngularMomentum(1,:);
AM2 = AngularMomentum(2,:);
AM3 = AngularMomentum(3,:);
%% Plotting Code for Part B
%{
figure(2);
subplot(2,1,1);
plot(t,[RR1,RR2,RR3,RR4,RR5,RR6,RR7,RR8,RR9])
title('R V. Time', 'FontSize', 18)
ylabel('Rotation Matrix')
xlabel('Time')
legend('R_1','R_2','R_3','R_4','R_5','R_6','R_7','R_8','R_9','Location','bestoutside')

subplot(2,1,2);
plot(t,[ww1,ww2,ww3])
title('Angular Velocity V. Time', 'FontSize', 18)
ylabel('Angular Velocity')
xlabel('Time')
legend('\omega_1','\omega_2','\omega_3','Location','bestoutside')

set(gcf,'paperorientation','portrait');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf',sprintf('plots.pdf',w0x));
%}
%% Plotting Code for Part C
%{
figure(2);
subplot(3,1,1);
plot(t,[RR1,RR2,RR3,RR4,RR5,RR6,RR7,RR8,RR9])
title('R V. Time', 'FontSize', 18)
ylabel('Rotation Matrix')
xlabel('Time')
legend('R_1','R_2','R_3','R_4','R_5','R_6','R_7','R_8','R_9','Location','bestoutside')

subplot(3,1,2);
plot(t,[ww1,ww2,ww3])
title('Angular Velocity V. Time', 'FontSize', 18)
ylabel('Angular Velocity')
xlabel('Time')
legend('\omega_1','\omega_2','\omega_3','Location','bestoutside')

subplot(3,1,3);
plot(t,[AM1',AM2',AM3'])
title('Angular Momentum V. Time', 'FontSize', 18)
ylabel('Angular Momentum')
xlabel('Time')
legend('h_1','h_2','h_3','Location','bestoutside')

set(gcf,'paperorientation','portrait');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf',sprintf('w0=1e%d.pdf',w0x));
%}
function xdot = f(t,x,J)
R = XtoR(x(1:9,1));
w = x(10:end,1);
zTorque = 100;
% %%
%
% YOUR CODE HERE TO COMPUTE Rdot and wdot
%
w1 = w(1);
w2 = w(2);
w3 = w(3);
what = [0 -w3 w2; w3 0 -w1; -w2 w1 0];

Rdot = R*what;

if t<7.578
    wd1 = ((J(2,2)-J(3,3))/(J(1,1)))*w3*w2;
    wd2 = ((J(3,3)-J(1,1))/(J(2,2)))*w3*w1;
    wd3 = (zTorque + (J(1,1)-J(2,2))*w1*w2)/J(3,3);
else
    wd1 = ((J(2,2)-J(3,3))/(J(1,1)))*w3*w2;
    wd2 = ((J(3,3)-J(1,1))/(J(2,2)))*w3*w1;
    wd3 = ((J(1,1)-J(2,2))*w1*w2)/J(3,3);
end
%Rdot = zeros(3);
wdot = [wd1;wd2;wd3];
%
% %%
xdot = [RtoX(Rdot); wdot];

function S = skew(w)
S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];

function X = RtoX(R)
X = reshape(R,9,1);

function R = XtoR(X)
R = reshape(X,3,3);

function [verts,faces,colors] = makebox(x,y,z)
verts = [0 x x 0 0 x x 0; 0 0 0 0 y y y y; 0 0 z z 0 0 z z];
verts = verts - repmat([x/2; y/2; z/2],1,size(verts,2));
faces = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 4 3 7 8; 5 6 2 1];
colors(:,:,1) = 1*ones(1,size(faces,1));
colors(:,:,2) = 1*ones(1,size(faces,1));
colors(:,:,3) = 0*ones(1,size(faces,1));
colors(1,5,1) = 1;
colors(1,5,2) = 0;
colors(1,5,3) = 0;
colors(1,3,1) = 0;
colors(1,3,2) = 0;
colors(1,3,3) = 1;
colors(1,2,1) = 0;
colors(1,2,2) = 1;
colors(1,2,3) = 0;



