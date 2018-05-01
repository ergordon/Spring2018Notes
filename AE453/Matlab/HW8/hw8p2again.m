function hw8p2again

Te = [0.5 0.6 0.7 0.8 0.9 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]; % Electron Temperature [eV]
OVi = [4.51e-25 3.02e-23 6.20e-22 6.04e-21 3.58e-20 1.5e-19 1.16e-17 1.08e-16 4.24e-16 1.08e-15 2.13e-15 3.59e-15 5.43e-15 7.61e-15 1.01e-14 1.28e-14 1.57e-14 1.88e-14 2.20e-14 2.53e-14 2.86e-14 3.20e-14 3.55e-14 3.90e-14];% Ionization Rate Coefficient [m^3/s]

pp = spline(Te,OVi);
qq = ppdiff(pp);
xx = linspace(Te(1),Te(end),200);
plot(xx,ppval(qq,xx))

end


function qq = ppdiff(pp,j)

if nargin < 1, help ppdiff, return, end
if nargin < 2, j = 1; end

% Check diff order
if ~isreal(j) || mod(j,1) || j < 0
    msgid = 'PPDIFF:DiffOrder';
    message = 'Order of derivative must be a non-negative integer!';
    error(msgid,message)
end

% Get coefficients
coefs = pp.coefs;
[m n] = size(coefs);

if j == 0
    % Do nothing
elseif j < n
    % Derivative of order J
    D = [n-j:-1:1; ones(j-1,n-j)];
    D = cumsum(D,1);
    D = prod(D,1);
    coefs = coefs(:,1:n-j);
    for k = 1:n-j
        coefs(:,k) = D(k)*coefs(:,k);
    end
else
    % Derivative kills PP
    coefs = zeros(m,1);
end

% Set output
qq = pp;
qq.coefs = coefs;
qq.order = size(coefs,2);
end