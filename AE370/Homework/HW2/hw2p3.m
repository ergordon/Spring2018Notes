function hw2p3
 t_0 = 0.0;
 t_f = 8.0;
 y_0 = [0.2, 0.0, 0.0].';
 h = 0.1;
 t = [t_0 : h : t_f];
 y_fe = forward_euler(t, y_0, @yprime,h);
 y_be = backward_euler(t, y_0, @yprime,h);
 y_fe = y_fe(1, :);
 y_be = y_be(1, :);
 y_a = 1.02 * exp(-4.0 * t) - 2.76667 * exp(-3.0 * t) ...
 + 2.00769 * exp(-2.0 * t) - 0.0610256 * cos(3.0 * t) ...
 - 0.0682051 * sin(3.0 * t);
 figure
 plot(t, y_a, 'r', t, y_fe, 'g', t, y_be, 'b', 'LineWidth', 2)
 set(gca, 'FontSize', 16)
 title(sprintf('Solutions using h = %4.2f ', h))
 xlabel('t')
 ylabel('y(t)')
 legend('Analytical', 'Forward Euler', 'Backward Euler')
 set(gcf,'paperorientation','landscape');
set(gcf,'paperunits','normalized');
set(gcf,'paperposition',[0 0 1 1]);
print(gcf,'-dpdf','h01.pdf');
end
function y = forward_euler(t, y_0, yp,h)
 %% Integrate the ODE given by yp using the forward Euler method.
 y = zeros(length(y_0), length(t)); %Create empty matrix to store all Y's
 y(:, 1) = y_0; %Intitial Condition
 for k = 1 : length(t) - 1
     y(:, k + 1) = y(:, k)+h*yp(t(k), y(:, k));
 end
end
function y = backward_euler(t, y_0, yp,h)
 %% Integrate the ODE given by yp using the backward Euler method.
 %% Fill in required function
 tol = 1.0e-6;
 y = zeros(length(y_0), length(t)); % Create empty matrix to store all Y's
 y(:, 1) = y_0;                     % Initial Condition
 for k = 1 : length(t) - 1
     y_k = y(:, k);                    % y_k   
     y_kp1_im1 = zeros(size(y_k));     %y_{k+1}^{i-1}
     y_kp1_i = y(:, k) + h * yp(t(k), y(:, k));   %y_{k+1}^{i}
     while norm(y_kp1_i - y_kp1_im1) > tol            % Test for Convergence
         y_kp1_im1 = y_kp1_i;                         %Update y_{k+1}^{i-1}
         y_kp1_i = y_k + h * yp(t(k + 1), y_kp1_im1); %Update y_{k+1}^{i}
     end
     y(:, k + 1) = y_kp1_i; 
     end
end
function yp = yprime(t, y)
 %% Evaluate the derivative y'(t, y).
 yp = zeros(3, 1);
 yp(1) = y(2);
 yp(2) = y(3);
 yp(3) = - 24.0 * y(1) - 26.0 * y(2) - 9.0 * y(3) + 7.0 * sin(3.0 * t);
end