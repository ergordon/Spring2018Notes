% Emilio Rey Gordon - Extended Kalman Filter
% User Guide
% 1. Choose Your Attitude Determination Instrument
%     a. 'Magnetometer' - Low Accuracy but Quick Update Time
%     b. 'EarthSensor' - Medium Accuracy, Quick Update Time
%     c. 'SunSensor' - Medium Accuracy, Slow Update Time
%     d. 'StarTracker' - High Accuracy, Slow Update Time
%     e. 'Perfect' - Experiences No Noise and Quick Update Time
% 2. Choose an end time.
% 3. Choose The Mode
%     a. 'NoControl' - Simply Track the state with no input (control) starting from initial conditions
%     b. 'ReferenceTracking' - Control is applied to establish a referenced input
% 4. Run ExtendedKalman('Instrument',tMax,'Mode')
% 
% The code will produce 3 figures:
%     1. Attitude V Time; Angular Velocity V. Time
%     2. Attitude Comparisons between Estimated state and observed state
%     3. Error V. Time
% 
% Example Operations
%   ExtendedKalman('StarTracker',100,'ReferenceTracking')
%   ExtendedKalman('Perfect',300,'NoControl')
%   ExtendedKalman('EarthSensor',50,'NoControl')
%
% Enjoy!

function ExtendedKalman(Instrument, tMax, Mode)

%Instrument Gallery: [Noise Error, Update Rate]
if strcmp(Instrument,'Magnetometer') == 1
    n = [1, 0.1]; % Magnetometer
elseif strcmp(Instrument,'EarthSensor') == 1
    n = [0.2, 0.1]; % Earth Sensor
elseif strcmp(Instrument,'SunSensor') == 1
	n = [0.15, 0.25]; % Sun Sensor
elseif strcmp(Instrument,'StarTracker') == 1
	n = [5.6e-4, 0.25]; % Star Tracker 
elseif strcmp(Instrument,'Perfect') == 1
	n = [0, 0.1]; % Perfect Instrument
else
    disp('Error')
end

wNoise = n(1);
hUpdateRate = n(2);

if strcmp(Mode,'NoControl') == 1
    [z, xhat, b, t] = NoControl(wNoise,hUpdateRate,tMax);
    plotter(z, xhat, b , t)
elseif strcmp(Mode,'ReferenceTracking') == 1
    [z, xhat, b, t] = ReferenceTracking(wNoise,hUpdateRate,tMax);
    plotter(z, xhat, b , t)
else
end
end

function [z, xhat, b, t] = NoControl(wNoise,hUpdateRate,tMax)
t = [0:hUpdateRate:tMax];

J1 = 12;
J2 = 14;
J3 = 8; 

syms w1 w2 w3 tt

%Observation Matrix
H = eye(6);

%State Matrix
F = [0 0 0        1               0               0        ;...
     0 0 0         0               1               0       ;...
     0 0 0         0               0               1       ;...
     0 0 0         0       ((J2-J3)/J1)*w3 ((J2-J3)/J1)*w2 ;...
     0 0 0 ((J3-J1)/J2)*w3         0       ((J3-J1)/J2)*w1 ;...
     0 0 0 ((J1-J2)/J3)*w2 ((J1-J2)/J3)*w1         0       ];

Q = eye(6);
R = hUpdateRate*eye(6);

% Initial Predictions
xhat(:,1) = [10; 0; 0; 1; 1; 1];    % Initial condition on the state
z(:,1) = H*xhat(:,1);               % Initial Observed State
P = Q;                              % Initial Error Covariance
b(:,1) = z-H*xhat(:,1);             % Initial Error

for i = 1:length(t)-1

  w = wNoise*(randn(1,1)-0.5); % Noise Generator

  % Update
  FF = vpa(subs(F,[w1 w2 w3],[xhat(4,i) xhat(5,i) xhat(6,i)]));
  
  K = P*H*inv(R); % Kalman Gain
  dP = FF*P + P*FF' + Q - K*R*K';
  dxhat = FF*xhat(:,i) + K*(z(:,i) - H*xhat(:,i));
  z(:,i+1) = H*xhat(:,i)+w; % Observed State

  xhat(:,i+1) = xhat(:,i) + hUpdateRate*dxhat; % State Estimate
  P = P + hUpdateRate*dP;
  b(:,i+1) = z(:,i)-H*xhat(:,i+1); %Error Update
end
end

function [z, xhat, b, t] = ReferenceTracking(wNoise,hUpdateRate,tMax)
t = [0:hUpdateRate:tMax];
  
J1 = 12; 
J2 = 14;
J3 = 8; 

syms w1 w2 w3 tt

%Observation Matrix
H = eye(6);

%State Matrix
F = [0 0 0        1               0               0        ;...
     0 0 0         0               1               0       ;...
     0 0 0         0               0               1       ;...
     -0.002*exp(0.009*tt) 0 0 0 0 0 ;...
     0 -0.002*exp(0.009*tt) 0 0 0 0 ;...
     0 0 -0.0002 0 0 0];
 
Q = eye(6);
R = eye(6);

% Initial Predictions
xhat(:,1) = [0; 10; 0; 0; 0; 0];    % Initial condition on the state
z(:,1) = H*xhat(:,1);               % Initial Observed State
P = Q;                              % Initial Error Covariance
b(:,1) = z-H*xhat(:,1);             % Initial Error

for i = 1:length(t)-1
  
  w = wNoise*(randn(1,1)-0.5); % Noise Generator
  ww(i) = w;
  
  % Update
  FF = vpa(subs(F,tt,t(i)));
  K = P*H*inv(R); % Kalman Gain
  dP = FF*P + P*FF' + Q - K*R*K';
  dxhat = FF*xhat(:,i) + K*(z(:,i) - H*xhat(:,i));
  z(:,i+1) = H*xhat(:,i)+w;

  xhat(:,i+1) = xhat(:,i) + hUpdateRate*dxhat; % State Estimate
  P = P + hUpdateRate*dP;
  b(:,i+1) = z(:,i)-H*xhat(:,i+1); %Error Update
end
end

function plotter(z,xhat,b,t)
%Estimated Values
theta1 = xhat(1,:);
theta2 = xhat(2,:);
theta3 = xhat(3,:);
omega1 = xhat(4,:);
omega2 = xhat(5,:);
omega3 = xhat(6,:);

%Observed Values
Otheta1 = z(1,:);
Otheta2 = z(2,:);
Otheta3 = z(3,:);
Oomega1 = z(4,:);
Oomega2 = z(5,:);
Oomega3 = z(6,:);

%Error
theta1Error = b(1,:);
theta2Error = b(2,:);
theta3Error = b(3,:);
omega1Error = b(4,:);
omega2Error = b(5,:);
omega3Error = b(6,:);

figure(1)
subplot(211)
plot(t,theta1)
hold on
plot(t,theta2)
plot(t,theta3)
title('Attitude V. Time')
legend('\theta_1','\theta_2','\theta_3')

subplot(212)
plot(t,omega1)
hold on
plot(t,omega2)
plot(t,omega3)
title('Angular Velocity V. Time')
legend('\omega_1','\omega_2','\omega_3')

figure(2)
subplot(211)
plot(t,theta1Error)
hold on
plot(t,theta2Error)
plot(t,theta3Error)
title('Attitude Error V. Time')
legend('\theta_1 Error','\theta_2 Error','\theta_3 Error')

subplot(212)
plot(t,omega1Error)
hold on
plot(t,omega2Error)
plot(t,omega3Error)
title('Angular Velocity Error V. Time')
legend('\omega_1 Error','\omega_2 Error','\omega_3 Error')

figure(3)
subplot(311)
plot(t,theta1)
hold on 
plot(t,Otheta1,'.')
title('\theta_1 V. Time')
legend('\theta_1 Estimated','\theta_1 Observed')

subplot(312)
plot(t,theta2)
hold on 
plot(t,Otheta2,'.')
title('\theta_2 V. Time')
legend('\theta_2 Estimated','\theta_2 Observed')

subplot(313)
plot(t,theta3)
hold on 
plot(t,Otheta3,'.')
title('\theta_3 V. Time')
legend('\theta_3 Estimated','\theta_3 Observed')

figure(4)
subplot(311)
plot(t,omega1)
hold on 
plot(t,Oomega1,'.')
title('\omega_1 V. Time')
legend('\omega_1 Estimated','\omega_1 Observed')

subplot(312)
plot(t,omega2)
hold on 
plot(t,Otheta2,'.')
title('\omega_2 V. Time')
legend('\omega_2 Estimated','\omega_2 Observed')

subplot(313)
plot(t,omega3)
hold on 
plot(t,Oomega3,'.')
title('\omega_3 V. Time')
legend('\omega_3 Estimated','\omega_3 Observed')
end













