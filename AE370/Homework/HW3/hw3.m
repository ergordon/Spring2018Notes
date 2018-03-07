function hw3
tol = 0.00001 %Tolerance for iterations stop
A = [5 -15 -1 2;...
     3 20 4 0;...
     1 6 11 -3;...
     4 11 3 12];
 
 B = [11; -8; 7; 4];
 
[uJ iterationsJ] = jacobiTolerance(A,B,tol);
[uGS iterationsGS] = gaussSeidelTolerance(A,B,tol);
exact = A\B;
endJ = uJ(:,end);
Jite = iterationsJ(end);
endGS = uGS(:,end);
GSite = iterationsGS(end);

sprintf('The Jacobi Method converges within a tolerance of %1.0d after the %d iteration with a final solution u=[%0.4f %0.4f %0.4f %0.4f]',tol,Jite,endJ(1),endJ(2),endJ(3),endJ(4))
sprintf('The Gauss-Seidel Method converges within a tolerance of %1.0d after the %d iteration with a final solution u=[%0.4f %0.4f %0.4f %0.4f]',tol,GSite,endGS(1),endGS(2),endGS(3),endGS(4))
sprintf('The Matlab Function A\\B = u has the final solution of u=[%0.4f %0.4f %0.4f %0.4f]',round(exact(1),5),round(exact(2),5),round(exact(3),5),round(exact(4),5))


end

function [uJ iterations] = jacobiTolerance(A,B,tolerance)
u = [0;0;0;0]; 
M = zeros(size(A));
N = -A;

for i=1:size(A,1)
    M(i,i) = A(i,i);
    N(i,i) = 0;
end
i = 0;
while norm(B-A*u)>tolerance
    i=i+1;
    u = inv(M)*(B+N*u);
    iterations(i) = i;
    uJ(1:4,i) = u;
end
end

function [uGS iterations] = gaussSeidelTolerance(A,B,tolerance)
u = [0;0;0;0]; 
M = triu(A)
N = M-A
i = 0;
%while norm(u-(A\B))>tolerance
while norm(B-A*u)>tolerance
    i=i+1;
    u = inv(M)*(B+N*u);
    iterations(i) = i;
    uGS(1:4,i) = u;
end
end

function plotter(A,B,tol)
[uJ iterationsJ] = jacobiTolerance(A,B,tol);
[uGS iterationsGS] = gaussSeidelTolerance(A,B,tol);

for j=1:length(uGS)
    GSRE(j) = norm(uGS(:,j)-(A\B))/norm(A\B);
    GSresidue(j) = norm(B-A*uGS(:,j));
end

for k=1:length(uJ)
    JRE(k) = norm(uJ(:,k)-(A\B))/norm(A\B);
    Jresidue(k) = norm(B-A*uJ(:,k));
end
figure(1)
plot(iterationsGS, GSRE,'LineWidth',2)
hold on
plot(iterationsJ, JRE,'LineWidth',2)
title('Error V. Iteration', 'FontSize', 18);
ylabel('Relative Error', 'FontSize', 13)
xlabel('Iteration', 'FontSize', 13)
grid on
grid minor
legend('Gauss-Seidel Method',...
    'Jacobi Method',...
    'Location','southeast')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','plots.pdf');

figure(2)
semilogy(iterationsGS,GSresidue,iterationsJ,Jresidue,'linewidth',2)
title('Jacobi V. Gauss-Seidel - Convergence', 'FontSize', 18);
ylabel('||B-Au^{(k)}||', 'FontSize', 13)
xlabel('Iteration Number (k)', 'FontSize', 13)
grid on
grid minor
legend('Gauss-Seidel Method',...
    'Jacobi Method',...
    'Location','southeast')

set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','plots2.pdf');
end
