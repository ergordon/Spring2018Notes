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
% Key parameters
L=1.0;    % x-dimension of domain (in m)
H=0.5;    % y-dimension of step (in m)
Nx= 2; %input(' Enter number of grid spacings in x (Nx) ')
Ny= 24; %input(' Enter number of grid spacings in y (Ny) ')
dx=L/Nx;  % grid spacing in x direction
dy=H/Ny;  % grid spacing in y direction
eta=dx/dy;  % grid spacing ratio
V=1;    % imposed inflow velocity (in m/s) (the outflow velocity = 2V)
rho=1;  % fluid density (in kg/m^3) (to compute the dynamic pressure)
%
% Compute global equation number of corners (see schematic above)
qA=1;
qB=2*Ny+1;
qC=Nx*(2*Ny+1)+1;
qD=qC+Ny;
qE=qD+Ny;
qF=qE+(Nx-1)*(Ny+1)+1;
qG=qF+Ny;
Numeq=qG     % number of equations
%
% Set up matrix (Amat) and vector (bvec) dimensions
Amat=sparse(Numeq,Numeq);
bvec=zeros(Numeq,1);
%
% Build linear system
for i=1:Numeq;       % place 1 along diagonal for all DOF then overwrite the interior grid points
    Amat(i,i)=1;
end
for i=2:Nx     % loop over interior grid points - left half of domain
    for j=2:2*Ny
        qij=compute_q(i,j,Nx,Ny);   % compute equation number q for (i,j) grid point
        Amat(qij,qij)=XXXXXXX;
        q=compute_q(i-1,j,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i+1,j,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i,j-1,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i,j+1,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
    end
end
for i=Nx+1:2*Nx     % loop over interior grid points - right half of domain
    for j=Ny+2:2*Ny
        qij=compute_q(i,j,Nx,Ny);
        Amat(qij,qij)= XXXXXXX;
        q=compute_q(i-1,j,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i+1,j,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i,j-1,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
        q=compute_q(i,j+1,Nx,Ny);
        Amat(qij,q)= XXXXXXX;
    end
end
%
% Build right-hand-side vector (imposed Psi BC)
for j=1:2*Ny+1     % loop over the left edge: 
    bvec(j)= XXXXXXX;
end
for i=2:2*Nx       % loop over top edge: 
    j=2*Ny+1;
    q=compute_q(i,j,Nx,Ny);
    bvec(q)= XXXXXXX;
end
for j=Ny+1:2*Ny+1  % loop over right edge: 
    i=2*Nx+1;
    q=compute_q(i,j,Nx,Ny);
    bvec(q)= XXXXXXX;
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
xvec=linspace(0,2*L,2*Nx+1);    % vector with 2*Nx+1 values of x
yvec=linspace(0,2*H,2*Ny+1);    % vector with 2*Ny+1 values of y
contourf(xvec,yvec,Psimat',15)  % create filled contour plots of phi field
xlabel('x (m)');
ylabel('y (m)');
title('Contour plot of function Phi');
%
% Compute and display velocity vector field
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
figure(2) % velocity vector plot
[x,y]=meshgrid(xvec,yvec);
quiver(x',y',vx,vy);    % create vector plot
xlabel('x (m)');
ylabel('y (m)');
title(' Velocity vector plot');
%
% Compute and display pressure field
figure(3)
pressure=0.5*rho*(vx.^2+vy.^2);    % pressure array
contourf(xvec,yvec,pressure',15)     % filled contour plot of pressure field
xlabel('x (m)');
ylabel('y (m)');
title(' Contour plot of pressure field');
%
% Display the pressure distribution along side CD (0<=y<=H)
figure(4)
yvalues=linspace(0,H,Ny+1);     % vector with Ny+1 values of y along CD
pvalues(1:Ny+1)=pressure(Nx+1,1:Ny+1);  % extract pressure values along CD from pressure array
plot(yvalues,pvalues,'bo-','linewidth',2)
xlabel('y (m)');
ylabel('pressure (Pa)');
title(' Pressure distribution along side CD');
end

function q=compute_q(i,j,Nx,Ny)
% Subroutine to compute the global equation number corresponding to grid
% point (i,j)
XXXXXXXXXXX
end
