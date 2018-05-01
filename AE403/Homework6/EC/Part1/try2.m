function try3
clc
m = 2;
l = 2;
g = 9.81;
h=.01;
tMax=2*pi;
%tMax=.2;
t = [0:h:tMax]';

A = [0 1 ; -g/l 0];
B = [0; 1];
C = [1 0];
D = 0;

%Define Gains Controller
Q = 1*[1 0; 0 1];
R = h;
K = lqr(A,B,Q,R)'


%% Initials
x = [pi/10; 0];     % Initial condition on the state
x0 = [pi/10; 0];    % Initial condition on the state
xhat0 = [pi/10; 0];
xhat = xhat0;
xreal=x0;
z(1) = C*xhat0;
P = Q;

for i = 1:length(t)-1
  w = 0.1*(randn(1,1)-0.5); % Noise
  
  xreal(:,i+1) = expm((A)*t(i))*x0;

  
  %Measured Value
  z(i+1) = C*xhat(:,i)+w;

  K = P*C'*inv(R);
  dxhat = A*xhat(:,i)+K*(z(i)-C*xhat(:,i)); %this is obersever. State Estimate
  dP = A*P+P*A'-P*C'*inv(R)*C*P+Q;
  

  P = P+h*dP;
  xhat(:,i+1) = xhat(:,i) + h*dxhat; % Next step for state estimate.

end
%{
figure(1)
subplot(211), plot(t,y,'--',t,ye,'-')
title('Time-varying Kalman filter response')
xlabel('Time'), ylabel('Output')
legend('y','yhat')
subplot(212), plot(t,y-z,'-.',t,y-ye,'-')
xlabel('Time'), ylabel('Output')
%}

thetaHat = xhat(1,:);
thetadotHat = xhat(2,:);
theta = xreal(1,:);
thetadot = xreal(2,:);

figure
plot(t,thetaHat,t,theta)
title('Time-varying Kalman filter response')
xlabel('Time'), ylabel('Theta')
legend('Estimated Theta','Actual Theta','location','southeast')
set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','EC1a.pdf')

figure
plot(t,thetaHat,t,theta,t,z,'.')
xlabel('Time'), ylabel('Theta')
legend('Estimated Theta','Actual Theta','Measured Theta w/ Noise','location','southeast')
set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','EC1b.pdf')
end