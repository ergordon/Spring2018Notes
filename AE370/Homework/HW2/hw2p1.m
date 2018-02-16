function hw2p1
% hw2p1 AE 370 Spring 2016 Homework 2 Problem 1
 x = linspace(0.0, 5.0, 101); % Plot f(x) and f'(x).
 f = fn(x);
 fp = fnp(x);
 figure
 plot(x, f, 'k', 'LineWidth', 2)
 grid on
 set(gca, 'FontSize', 18)
 title('Given function f(x)')
 xlabel('x')
 ylabel('f(x)')
 print(gcf, '-depsc2', 'hw2p1-fn')

 x_low = 1.0;
 x_high = 4.0;
 x_tol = 1.0e-6;
 x_hat_bis = bisection(@fn, x_low, x_high, x_tol);
 x_0 = x_low + 0.5 * (x_high - x_low); % Starting guess
 f_tol = 1.0e-6;
 x_hat_nr = newton_raphson(@fn, @fnp, x_0, f_tol);
 x_1 = x_0 - 0.01 * (x_high - x_low); % Second starting guess
 x_hat_sec = secant(@fn, x_0, x_1, f_tol);

 fprintf('Bisection: x_hat = %10.6f\n', x_hat_bis)
 fprintf('Newton-Raphson: x_hat = %10.6f\n', x_hat_nr)
 fprintf('Secant method: x_hat = %10.6f\n', x_hat_sec)
end
function f = fn(x)
% Evaluate the given function.
 f = 2.0 * x.^3 + 5.875 * x.^2 - 8.625 * x - 24.75;
end
function fp = fnp(x)
% Evaluate the derivative of the given function.
 fp = 6.0 * x.^2 + 11.75 * x - 8.625;
end
function x = bisection(f, x_low, x_high, x_tol)
% Isolate a root of f(x) using bisection.
% fill up this code here
while x_high-x_low > x_tol
    c = (x_low+x_high)/2;
    if fn(x_low)*fn(c) > 0
        x_low = c;
    else
        x_high = c;
    end
end
x = c;
end
function x = newton_raphson(f, fp, x_0, tol)
% Isolate a root of f(x) using Newton-Raphson iteration.
% fill up this code here
x = x_0;
while abs(fn(x))>tol
    dx = -fn(x)/fnp(x);
    x = x+dx;
end
end
function x = secant(f, x_0, x_1, tol)
% Isolate a root of f(x) using the secant method.
% fill up this code here
x = x_0;
while abs(fn(x))>tol
    xx = x_1 - fn(x_1)*((x_1 - x)/(fn(x_1)-fn(x)));
    x = x_1;
    x_1 = xx;
end
end