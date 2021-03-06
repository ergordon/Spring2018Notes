function hw2_prob8

% DEFINE THE TIME INTERVAL
tMax = 6;
tRate = 30;
t = linspace(0,tMax,tMax*tRate);

%DEFINE THE INITIAL CONDITIONS
%xH0 = [0;pi/2;0];
%R0 = [0 0 1; 0 1 0; -1 0 0];
%xQ0 = [1/sqrt(2);0;1/sqrt(2);0];
xH0 = [0;0;0];
R0 = eye(3);
xQ0 = [0; 0; 0; 1];
%Could use Identity matrix
QtoR(xQ0(1),xQ0(2:end))
HtoR(xH0)

% INTEGRATE THE ANGULAR RATE EQUATIONS
[t,xR] = ode45(@fR,t,RtoX(R0));
[t,xH] = ode45(@fH,t,xH0);
[t,xQ] = ode45(@fQ,t,xQ0);

%finalR = xR(end,:)
finalR = reshape(xR(end,:),3,3)
finalH = reshape(xH(end,:),3,1)
finalQqo = xQ(end,1)
finalQq = reshape(xQ(end,2:4),3,1)

finalR = reshape(xR(end,:),3,3)
reshape(QtoR(finalQqo,finalQq),3,3)
reshape(HtoR(finalH),3,3)

% PLOT THE States
figure(3);
clf;
subplot(3,1,1);
plot(t,xR);
title('R V. Time', 'FontSize', 18)
ylabel('Rotation Matrix')
xlabel('Time')
legend('R_1','R_2','R_3','R_4','R_5','R_6','R_7','R_8','R_9','Location','bestoutside')
subplot(3,1,2);
plot(t,xH);
title('Theta V. Time', 'FontSize', 18)
ylabel('Euler Angle')
xlabel('Time')
legend('\theta_1','\theta_2','\theta_3','Location','bestoutside')
subplot(3,1,3);
plot(t,xQ);
title('Q V. Time', 'FontSize', 18)
ylabel('Quaternion')
xlabel('Time')
legend('q_0','q_1','q_2','q_3','Location','bestoutside')
set(gcf,'paperorientation','portrait');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','figure3.pdf');


% FIND THE RATES
for i=1:length(t)
    rdot(:,i) = RtoX(getRdot(t(i),XtoR(xR(i,:)')));
    thetadot(:,i) = getHdot(t(i),xH(i,:)');
    [q0dot(i),qdot(:,i)] = getQdot(t(i),xQ(i,1),xQ(i,2:4)');
end

% PLOT THE RATES
figure(1);
clf;
subplot(3,1,1);
plot(t,rdot);
title('Rdot V. Time', 'FontSize', 18)
ylabel('d/dt Rotation Matrix')
xlabel('Time')
legend('dR_1','dR_2','dR_3','dR_4','dR_5','dR_6','dR_7','dR_8','dR_9','Location','bestoutside')
subplot(3,1,2);
plot(t,thetadot);
title('Thetadot V. Time', 'FontSize', 18)
ylabel('d/dt Euler Angle')
xlabel('Time')
legend('d\theta_1','d\theta_2','d\theta_3','Location','bestoutside')
subplot(3,1,3);
plot(t,q0dot,t,qdot);
title('Qdot V. Time', 'FontSize', 18)
ylabel('d/dt Quaternion')
xlabel('Time')
legend('dq_0','dq_1','dq_2','dq_3','Location','bestoutside')
set(gcf,'paperorientation','portrait');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','figure1.pdf');

% CREATE A BOX
[p1,faces,colors] = makebox(0.5,1,0.2);

% SETUP THE PLOT 
figure(2);
clf;
axis(1.25*[-1 1 -1 1 -1 1]);
%axis equal;
hold on;
plot3(0,0,0,'k.','markersize',16);
hBox = patch('Vertices',p1','Faces',faces,...
          'CData',colors,'FaceColor','flat');
hTitle = title(sprintf('t = %4.2f',0));
lighting flat
light('Position',[0 -2 -1])
light('Position',[0 -2 1])

% ANIMATE THE RESULTS
x = xR;
firsttime = 1;
i = 1;
dt = max(t)/(length(t)-1);
tic;
while (i<length(t))
    if (toc > dt)
        tic;
        i = i+1;
        R = XtoR(x(i,:));
        p0 = R*p1;
        set(hBox,'Vertices',p0');
        set(hTitle,'string',sprintf('t = %4.2f',t(i)));
        
        if (firsttime)
            firsttime=0;
            pause(0.5);
        end
        drawnow;
    end
end

function R=HtoR(H)
% CODE TO CONVERT FROM XYZ EULER ANGLES TO A ROTATION MATRIX
theta1 = H(1);
theta2 = H(2);
theta3 = H(3);

R11 = cos(theta2)*cos(theta3);
R12 = -cos(theta2)*sin(theta3);
R13 = sin(theta2);

R21 = sin(theta1)*sin(theta2)*cos(theta3)+cos(theta1)*sin(theta3);
R22 = -sin(theta1)*sin(theta2)*sin(theta3)+cos(theta1)*cos(theta3);
R23 = -sin(theta1)*cos(theta2);

R31 = -cos(theta1)*sin(theta2)*cos(theta3)+sin(theta1)*sin(theta3);
R32 = cos(theta1)*sin(theta2)*sin(theta3)+sin(theta1)*cos(theta3);
R33 = cos(theta1)*cos(theta2);

R=[R11 R12 R13; R21 R22 R23; R31 R32 R33];


function R=QtoR(q0,q)
% CODE TO CONVERT FROM A QUATERNION TO A ROTATION MATRIX
qhat = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0];

R = ((q0^2)-q'*q)*eye(3)+2*q*q'+2*q0*qhat;


function [q0dot, qdot] = getQdot(t,q0,q)
% CODE TO COMPUTE Qdot = (q0dot,qdot)
omega = 10*exp(-t)*[sin(t);sin(2*t);sin(3*t)];
omegahat = [0 -omega(3) omega(2); omega(3) 0 -omega(1); -omega(2) omega(1) 0];
       
q0dot = -0.5*omega'*q;
qdot = 0.5*((omega*q0)-(omegahat*q));


function thetadot = getHdot(t,theta)
% CODE TO COMPUTE Thetadot
theta2 = theta(2);
theta3 = theta(3);
omega = 10*exp(-t)*[sin(t);sin(2*t);sin(3*t)];

B=[(cos(theta3)/cos(theta2)) (-sin(theta3)/cos(theta2)) 0;...
    (sin(theta3)) (cos(theta3)) 0;...
    (-tan(theta2)*cos(theta3)) (tan(theta2)*sin(theta3)) 1];

thetadot = B*omega;

function Rdot = getRdot(t,R)
% CODE TO COMPUTE Rdot
omega = 10*exp(-t)*[sin(t);sin(2*t);sin(3*t)];
omegahat = [0 -omega(3) omega(2); omega(3) 0 -omega(1); -omega(2) omega(1) 0];

Rdot = R*omegahat;


function xdot = fQ(t,x)
q0 = x(1,1);
q = x(2:4,1);
[q0dot, qdot] = getQdot(t,q0,q);
xdot = [q0dot; qdot];

function xdot = fH(t,x)
xdot = getHdot(t,x);

function xdot = fR(t,x)
R = XtoR(x);
xdot = RtoX(getRdot(t,R));

function X = RtoX(R)
X = reshape(R,9,1);

function R = XtoR(X)
R = reshape(X,3,3);

function [verts,faces,colors] = makebox(x,y,z)
verts = [0 x x 0 0 x x 0; 0 0 0 0 y y y y; 0 0 z z 0 0 z z];
faces = [1 2 3 4; 2 6 7 3; 6 5 8 7; 5 1 4 8; 4 3 7 8; 5 6 2 1];
colors(:,:,1) = 1*ones(1,size(faces,1));
colors(:,:,2) = 1*ones(1,size(faces,1));
colors(:,:,3) = 0*ones(1,size(faces,1));
colors(1,5,1) = 1;
colors(1,5,2) = 0;
colors(1,5,3) = 0;