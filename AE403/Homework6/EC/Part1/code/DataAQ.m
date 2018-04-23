% Loop over each flight

DesignProblem03('Controller','datafile','data.mat','display',false);
% Load data
load('data.mat')

% Get t and x
t = processdata.t;
x = processdata.x;

xdotHat
zdotHat 
thetadotHat
thetaHat
phiHat

plot(tt, x)
hold on

%{
    set(gca,'fontsize',14);
    xlabel('Distance (x)');
    ylabel('Frequency');
    set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    title('Frequency Distribution of 1000 Simulated Flights')
    print(gcf,'-dpdf','flights.pdf');
mean(Xf)
median(Xf)

m = matfile(xFlight1,'Writable',isWritable)
save(xFlight,Xf)
%}