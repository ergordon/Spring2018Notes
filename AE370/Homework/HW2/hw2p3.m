function hw2p3
 t_0 = 0.0;
 t_f = 8.0;
 y_0 = [0.2, 0.0, 0.0].';
 h = 0.01;
 t = [t_0 : h : t_f];
 y_fe = forward_euler(t, y_0, @yprime);
 y_be = backward_euler(t, y_0, @yprime);
 y_fe = y_fe(1, :);
 y_be = y_be(1, :);
 y_a = 1.02 * exp(-4.0 * t) - 2.76667 * exp(-3.0 * t) ...
 + 2.00769 * exp(-2.0 * t) - 0.0610256 * cos(3.0 * t) ...
 - 0.0682051 * sin(3.0 * t);
 figure
 plot(t, y_a, 'r', t, y_fe, 'g', t, y_be, 'b', 'LineWidth', 2)
 set(gca, 'FontSize', 18)
 title(sprintf('Solutions using h = %4.2f ', h))
 xlabel('t')
 ylabel('y(t)')
 legend('Analytical', 'Forward Euler', 'Backward Euler')
 print(gcf, '-depsc2', 'hw5p1-solns')
end
function y = forward_euler(t, y_0, yp)
 %% Integrate the ODE given by yp using the forward Euler method.
 %% Fill in required function
end
function y = backward_euler(t, y_0, yp)
 %% Integrate the ODE given by yp using the backward Euler method.
 %% Fill in required function
end
function yp = yprime(t, y)
 %% Evaluate the derivative y'(t, y).
 yp = zeros(3, 1);
 yp(1) = y(2);
 yp(2) = y(3);
 yp(3) = - 24.0 * y(1) - 26.0 * y(2) - 9.0 * y(3) + 7.0 * sin(3.0 * t);
end