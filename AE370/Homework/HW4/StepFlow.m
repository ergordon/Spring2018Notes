% Program StepFlow
% Finite difference code to solve the problem of an inviscid,
% incompressible flow going over a step.
% 
%
%     B                        E                       G
%     --------------------------------------------------
%     |                                                |
%     |                                                |
%     |                                                |
%     |                                                |
%     |                                                |
%     |                         D                      |F
%     |                         ------------------------
%     |                         |
%     |                         |
%     |                         |
%     |                         |
%     |                         |
%     |                         |
%     ---------------------------
%    A                          C
%
% Dimensions:    AC = DF = BE = EG = L
%                CD = FG = DE = AB/2 = H
% The domain is discretized with Nx grid spacings along AC and DF
% and NY grid spacings along CD and FG.
% Index convention:
%     Local indices (i,j)         Global index (q)
%   A      (1,1)                         qA=1
%   B      (1,2*Ny+1)                    qB=2*Ny+1
%   C      (Nx+1,1)                      qC=Nx*(2*Ny+1)+1
%   D      (Nx+1,Ny+1)                   qD=qC+Ny
%   E      (Nx+1,2*Ny+1)                 qE=qD+Ny
%   F      (2*Nx+1,Ny+1)                 qF=qE+(Nx-1)*(Ny+1)+1
%   G      (2*Nx+1,2*Ny+1)               qG=qF+Ny
%
function StepFlow
close all; clear all; clc;
% Key parameters
L = 1; %x-dimension of domain (m)
H = 0.2; %y-direction of step (m)

Ny1= 20; %input(' Enter number of grid spacings in y (Ny) ')
Nx1= 2*Ny1; %input(' Enter number of grid spacings in x (Nx) ')
Ny2 = 40;
Nx2 = 2*Ny2;
Ny3 = 80;
Nx3 = 2*Ny3;
Ny4 = 160;
Nx4 = 2*Ny4;

dx1=L/Nx1;  % grid spacing in x direction
dy1=H/Ny1;  % grid spacing in y direction
dx2=L/Nx2;  
dy2=H/Ny2;  
dx3=L/Nx3;  
dy3=H/Ny3;  
dx4=L/Nx4;  
dy4=H/Ny4;  

V=1;    % imposed inflow velocity (in m/s) (the outflow velocity = 2V)
rho=1;  % fluid density (in kg/m^3) (to compute the dynamic pressure)

Psimat1 = computeStream(Nx1, Ny1,L,H,V,dx1,dy1);
Psimat2 = computeStream(Nx2, Ny2,L,H,V,dx2,dy2);
Psimat3 = computeStream(Nx3, Ny3,L,H,V,dx3,dy3);
Psimat4 = computeStream(Nx4, Ny4,L,H,V,dx4,dy4);

plotContour(Psimat4,Nx4,Ny4,L,H)

figure(2)
plotVectorField(Psimat1,Nx1,Ny1,L,H,V,dx1,dy1)

figure(3)
subplot(2, 2,1);
plotContourPressures(Psimat1,Nx1,Ny1,L,H,V,dx1,dy1,rho)
subplot(2, 2,2);
plotContourPressures(Psimat2,Nx2,Ny2,L,H,V,dx2,dy2,rho)
subplot(2, 2,3);
plotContourPressures(Psimat3,Nx3,Ny3,L,H,V,dx3,dy3,rho)
subplot(2, 2,4);
plotContourPressures(Psimat4,Nx4,Ny4,L,H,V,dx4,dy4,rho)

figure(5)
subplot(2, 2,1);
plotContourPressuresZoomed(Psimat1,Nx1,Ny1,L,H,V,dx1,dy1,rho)
subplot(2, 2,2);
plotContourPressuresZoomed(Psimat2,Nx2,Ny2,L,H,V,dx2,dy2,rho)
subplot(2, 2,3);
plotContourPressuresZoomed(Psimat3,Nx3,Ny3,L,H,V,dx3,dy3,rho)
subplot(2, 2,4);
plotContourPressuresZoomed(Psimat4,Nx4,Ny4,L,H,V,dx4,dy4,rho)

figure(4)
plotPressures(Psimat1,Psimat2, Psimat3, Psimat4,Nx1,Ny1,Nx2,Ny2,Nx3,Ny3,Nx4,Ny4,H,V,dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4,rho)

end

function [Psimat] = computeStream(Nx, Ny, L, H,V,dx,dy)
eta=dx/dy;  % grid spacing ratio
rho=1;  % fluid density (in kg/m^3) (to compute the dynamic pressure)

% Compute global equation number of corners (see schematic above)
qA=1;
qB=2*Ny+1;        % N
qC=Nx*(2*Ny+1)+1; % Nx*(qB)+1
qD=qC+Ny;
qE=qD+Ny;
qF=qE+(Nx-1)*(Ny+1)+1;
qG=qF+Ny;
Numeq=qG;     % number of equations
%
% Set up matrix (Amat) and vector (bvec) dimensions
Amat=sparse(Numeq,Numeq);
bvec=zeros(Numeq,1);
% Build linear system
for i=1:Numeq;       % place 1 along diagonal for all DOF then overwrite the interior grid points
    Amat(i,i)=1;
end
for i=2:Nx % loop over interior grid points - left half of domain
    for j=2:2*Ny
        qij=compute_q(i,j,Nx,Ny); % compute equation number q for(i,j) grid point 
        Amat(qij,qij)=-2*(1+eta^2);%Grid Point
        q=compute_q(i-1,j,Nx,Ny); %Grid Point to the left
        Amat(qij,q)= 1;
        q=compute_q(i+1,j,Nx,Ny); %Grid Point to the right
        Amat(qij,q)= 1;
        q=compute_q(i,j-1,Nx,Ny); % Grid point below
        Amat(qij,q)= eta^2;
        q=compute_q(i,j+1,Nx,Ny); %Grid Point above
        Amat(qij,q)= eta^2;
    end
end
for i=Nx+1:2*Nx % loop over interior grid points - right half of domain
    for j=Ny+2:2*Ny
        qij=compute_q(i,j,Nx,Ny);
        Amat(qij,qij)= -2*(1+eta^2);
        q=compute_q(i-1,j,Nx,Ny);
        Amat(qij,q)= 1;
        q=compute_q(i+1,j,Nx,Ny);
        Amat(qij,q)= 1;
        q=compute_q(i,j-1,Nx,Ny);
        Amat(qij,q)= eta^2;
        q=compute_q(i,j+1,Nx,Ny);
        Amat(qij,q)= eta^2;
    end
end


%
% Build right-hand-side vector (imposed Psi BC)
for j=1:2*Ny+1 % loop over the left edge:
    bvec(j)= V*(j-1)*dy;
end
for i=2:2*Nx % loop over top edge:
    j=2*Ny+1;
    q=compute_q(i,j,Nx,Ny);
    bvec(q)= 2*V*H;
end
for j=Ny+1:2*Ny+1 % loop over right edge:
    i=2*Nx+1;
    q=compute_q(i,j,Nx,Ny);
    y = (j-1)*dy;
    bvec(q)= 2*V*(y-H);
end
% Remainder of boundary has Psi=0
%
% Solve linear system
Psivec=Amat\bvec;
% Build Psi array for visualization (fill lower right with zero
Psimat=zeros(2*Nx+1,2*Ny+1);
for i=1:2*Nx+1     % loop over all points in domain
    for j=1:2*Ny+1
        q=compute_q(i,j,Nx,Ny);
        if q == 0    % if inside step region, assign Phi=0
            Psimat(i,j)=0;
        else  % if inside computational domain, extra Phi value from solution vector
            Psimat(i,j)=Psivec(q);
        end
    end
end
end

function q=compute_q(i,j,Nx,Ny)
% Subroutine to compute the global equation number corresponding to grid
% point (i,j)
q = 0;
if(i <= Nx +1)
    q = (i-1)*(2*Ny+1)+j;
elseif(j >= Ny + 1)
    q = (Nx + 1)*(2*Ny+1)+(i-Nx-2)*(Ny+1)+j-Ny;
end
end

function plotContour(Psimat,Nx,Ny,L,H)

xvec=linspace(0,2*L,2*Nx+1);    % vector with 2*Nx+1 values of x
yvec=linspace(0,2*H,2*Ny+1);    % vector with 2*Ny+1 values of y
contourf(xvec,yvec,Psimat',15)  % create filled contour plots of phi field
colormap(jet);
xlabel('x (m)');
ylabel('y (m)');
title('Contour plot of function Phi');
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','contourStream.pdf');
end

function [vx,vy] = Velocity(Psimat,Nx,Ny,V,dx,dy)
vx=zeros(2*Nx+1,2*Ny+1);    % x-component of velocity at each grid point
vy=zeros(2*Nx+1,2*Ny+1);    % y-component of velocity at each grid point

for j=1:2*Ny+1     % left edge
    vx(1,j)=V;
    vy(1,j)=0;
end
for j=Ny+1:2*Ny+1   % right edge
    vx(2*Nx+1,j)=2*V;
    vy(2*Nx+1,j)=0;
end
for i=2:2*Nx       % top edge (backward difference)
    vx(i,2*Ny+1)=(Psimat(i,2*Ny+1)-Psimat(i,2*Ny))/dy;
    vy(i,2*Ny+1)=0;
end
for i=2*Nx         % bottom edge (forward difference)
    vx(i,1)=(Psimat(i,2)-Psimat(i,1))/dy;
end
for i=2:2*Nx         % interior nodes
    for j=2:2*Ny
        vx(i,j)=(Psimat(i,j+1)-Psimat(i,j-1))/2/dy;
        vy(i,j)=-(Psimat(i+1,j)-Psimat(i-1,j))/2/dx;
    end
end
end


function plotVectorField(Psimat,Nx,Ny,L,H,V,dx,dy)
% Compute and display velocity vector field

xvec=linspace(0,2*L,2*Nx+1);    % vector with 2*Nx+1 values of x
yvec=linspace(0,2*H,2*Ny+1);

[vx,vy]=Velocity(Psimat,Nx,Ny,V,dx,dy);

[x,y]=meshgrid(xvec,yvec);
quiver(x',y',vx,vy);    % create vector plot
xlabel('x (m)');
ylabel('y (m)');
axis([0 2 0 1])
title(' Velocity vector plot');
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','vfield.pdf');
end


function plotContourPressures(Psimat,Nx,Ny,L,H,V,dx,dy,rho)
[vx,vy] = Velocity(Psimat,Nx,Ny,V,dx,dy);
% Compute and display pressure field
xvec=linspace(0,2*L,2*Nx+1);    % vector with 2*Nx+1 values of x
yvec=linspace(0,2*H,2*Ny+1);

pressure=0.5*rho*(vx.^2+vy.^2);    % pressure array
contourf(xvec,yvec,pressure',30)     % filled contour plot of pressure field
colormap(jet);
colorbar
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Pressure Field w/ Ny = %d',Ny));
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','contourPressure.pdf');
end

function plotContourPressuresZoomed(Psimat,Nx,Ny,L,H,V,dx,dy,rho)
[vx,vy] = Velocity(Psimat,Nx,Ny,V,dx,dy);
% Compute and display pressure field
xvec=linspace(0,2*L,2*Nx+1);    % vector with 2*Nx+1 values of x
yvec=linspace(0,2*H,2*Ny+1);

for i = 1:length(xvec)
    if xvec(i) == 0.8
        xLow = i
    end
    if xvec(i) == 1.3
        xHigh = i
    end
end
for i = 1:length(yvec)
    if yvec(i) == 0.4
        yLow = i
    end
    if yvec(i) == 0.6
        yHigh = i
    end
end
pressure=0.5*rho*(vx.^2+vy.^2);    % pressure array
contourf(xvec(xLow:xHigh),yvec(yLow:yHigh),pressure(xLow:xHigh,yLow:yHigh)',(Ny/2)+5)     % filled contour plot of pressure field
colormap(jet);
colorbar
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Pressure Field w/ Ny = %d',Ny));
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','contourPressureZoomed.pdf');
end

function plotPressures(Psimat1,Psimat2, Psimat3, Psimat4,Nx1,Ny1,Nx2,Ny2,Nx3,Ny3,Nx4,Ny4,H,V,dx1,dy1,dx2,dy2,dx3,dy3,dx4,dy4,rho)
[vx1,vy1] = Velocity(Psimat1,Nx1,Ny1,V,dx1,dy1);
[vx2,vy2] = Velocity(Psimat2,Nx2,Ny2,V,dx2,dy2);
[vx3,vy3] = Velocity(Psimat3,Nx3,Ny3,V,dx3,dy3);
[vx4,vy4] = Velocity(Psimat4,Nx4,Ny4,V,dx4,dy4);

pressure1=0.5*rho*(vx1.^2+vy1.^2)
pressure2=0.5*rho*(vx2.^2+vy2.^2);
pressure3=0.5*rho*(vx3.^2+vy3.^2);
pressure4=0.5*rho*(vx4.^2+vy4.^2);

yvalues1=linspace(0,H,Ny1+1);     % vector with Ny+1 values of y along CD
yvalues2=linspace(0,H,Ny2+1);
yvalues3=linspace(0,H,Ny3+1);
yvalues4=linspace(0,H,Ny4+1);

pvalues1(1:Ny1+1)=pressure1(Nx1+1,1:Ny1+1)  % extract pressure values along CD from pressure array
pvalues2(1:Ny2+1)=pressure2(Nx2+1,1:Ny2+1);
pvalues3(1:Ny3+1)=pressure3(Nx3+1,1:Ny3+1);
pvalues4(1:Ny4+1)=pressure4(Nx4+1,1:Ny4+1);

plot(yvalues1,pvalues1,'bo-','linewidth',2)
hold on
plot(yvalues2,pvalues2,'rs-','linewidth',2)
hold on
plot(yvalues3,pvalues3,'g--','linewidth',2)
hold on
plot(yvalues4,pvalues4,'m+-','linewidth',2)
xlabel('y (m)');
ylabel('pressure (Pa)');
%axis([0.19 0.2 0 6])
title(' Pressure distribution along side CD');
legend(sprintf('N_y = %d',Ny1),sprintf('N_y = %d',Ny2),sprintf('N_y = %d',Ny3),sprintf('N_y = %d',Ny4),'Location','northwest')
set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','PressurePlots.pdf');
end