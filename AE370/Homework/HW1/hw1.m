%{
clc; clear;
X = [0; sqrt(3/5); -sqrt(3/5)];

Y = [0; sqrt(3/5); -sqrt(3/5)];

A11 = X(1)*sin((X(1))^2)+log(2+Y(1));
A12 = X(1)*sin((X(1))^2)+log(2+Y(2));
A13 = X(1)*sin((X(1))^2)+log(2+Y(3));

A21 = X(2)*sin((X(2))^2)+log(2+Y(1));
A22 = X(2)*sin((X(2))^2)+log(2+Y(2));
A23 = X(2)*sin((X(2))^2)+log(2+Y(3));

A31 = X(3)*sin((X(3))^2)+log(2+Y(1));
A32 = X(3)*sin((X(3))^2)+log(2+Y(2));
A33 = X(3)*sin((X(3))^2)+log(2+Y(3));

W = [(64/81) (40/81) (40/81); (40/81) (25/81) (25/81); (40/81) (25/81) (25/81)];

I33 = W(1,1)*A11+ W(1,2)*A12+ W(1,3)*A13+ W(2,1)*A21+ W(2,2)*A22+ W(2,3)*A23+ W(3,1)*A31+ W(3,2)*A32+ W(3,3)*A33;

I11 = 2.77258872;
I22 = 2.59858;

Iex = 2.59167;

(1-(I33/Iex))*100
%}
clc; clear;

N = 1:30;
intLB = 0; % Integral Lower Bound
intUB = 2*pi; %Integral Upper Bound

%syms x real
%f(x) = (cos(x^3 + 3*x))/(exp(x^2));

IRECT = [];
ITRAP = [];
for j=1:length(N)
    n = N(j);
    h = (intUB - intLB)/n; %Length of each SubInterval
    xi = intLB; % First x position
    for i=1:n
        xi = [xi xi(i)+h];
    end
    %% Rectangule Rule
    Irect = 0;
    for i=1:n
        xx = (xi(i)+xi(i+1))/2
        fxx = (cos(xx^3 + 3*xx))/(exp(xx^2));
        Irect = Irect + h*fxx;
    end
    IRECT = [IRECT Irect]
    %% Trapezoid Rule
    Itrap = 0;
    for i=1:n
        fx1 =(cos(xi(i)^3 + 3*xi(i)))/(exp(xi(i)^2));
        fx2 = (cos(xi(i+1)^3 + 3*xi(i+1)))/(exp(xi(i+1)^2));
        Itrap = Itrap + (h/2)*(fx1+fx2);
    end
    ITRAP = [ITRAP Itrap];
    j
end
plot(N,IRECT,('.-'),N,ITRAP,('x-'))
     p(1).LineWidth = 2;
     p(2).LineWidth = 2;
     title('Integration Method V. Number of Intervals', 'FontSize', 24)
     ylabel('Integral Value', 'FontSize', 18)
     xlabel('Number of Intervals', 'FontSize', 18)
     grid on
     grid minor
     xlim([0,30]);
     ylim([-1,3.5]);
     legend('Rectangle','Trapezoide')
     set(gcf,'paperorientation','landscape');
     set(gcf,'paperunits','normalized');
     set(gcf,'paperposition',[0 0 1 1]);
     print(gcf,'-dpdf','graph2.pdf');