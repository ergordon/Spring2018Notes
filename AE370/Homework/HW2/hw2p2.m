function hw2p2
clc;
x = [1, 1].';
x_hist = [x]; % Solution history, stored by cols
f_tol = 1.0e-6;
f = fn(x);
while norm(f) > f_tol
    x
    inv(jac(x))
    f
    Delta_x_other = - inv(jac(x))*f
    Delta_x = - jac(x) \ f
    x = x + Delta_x
    f = fn(x)
    x_hist(:, end+1) = x
end
fprintf('Found root x_1 = %10.6f, x_2 = %10.6f\n', x(1), x(2))
    function f = fn(x)
    % Evalute the given function.
    f = zeros(2, 1);
    f(1) = x(1)^2 + x(2) - 4.5;
    f(2) = sqrt(3)*x(1) -x(2)^2 - 0.75;
    end
    
    function jac = jac(x)
    % Evaluate the Jacobian of the given function.
    jac = zeros(2, 2);
    jac(1, 1) = 2.0 * x(1); 
    jac(1, 2) = 1;
    jac(2, 1) = sqrt(3); 
    jac(2, 2) = -2*x(2);
    end
end
