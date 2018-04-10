clc; clear;

A = -1:.001:1;

% get 2-D mesh for x and y
[a1 a3] = meshgrid(A);	

% check conditions for these values
cond1 = (1+(a1.*3)+(a1.*a3.*1)) > 0;
cond2 = ((1+(a1.*3)+(a1.*a3.*1)).^2 - 16*a1.*a3.*1) > 0;
cond3 = a1.*a3 > 0;
cond4 = a1 > a3;
cond5 = abs(a1) < 1;
cond6 = abs(a3) < 1;

% convert to double for plotting
cond1 = double(cond1);
cond2 = double(cond2);
cond3 = double(cond3);
cond4 = double(cond4);
cond5 = double(cond5);
cond6 = double(cond6);

% set the 0s to NaN so they are not plotted
cond1(cond1 == 0) = NaN;
cond2(cond2 == 0) = NaN;
cond3(cond3 == 0) = NaN;
cond4(cond4 == 0) = NaN;
cond5(cond5 == 0) = NaN;
cond6(cond6 == 0) = NaN;

% multiply the condaces to keep only the common points
cond = cond1.*cond2.*cond3.*cond4.*cond5.*cond6;

s = surf(a1,a3,cond);
s.EdgeColor = 'none';
s.LineStyle = '-';
s.FaceColor = 'flat';
s.FaceLighting = 'flat';
axis on;
axis([-1 1 -1 1])
    xlabel('a_1','FontSize',14)
    ylabel('a_3','FontSize',14)
    title(sprintf('Criteria of J1, J2, J3 to Guarantee (Marginal) Passive Gravity-Gradient Stabilization \n For a Spacecraft in Low-Earth Orbit.'))

set(gcf,'paperorientation','landscape');
    set(gcf,'paperunits','normalized');
    set(gcf,'paperposition',[0 0 1 1]);
    print(gcf,'-dpdf','ExtraCredit.pdf')
view(0,90)% change to top view
