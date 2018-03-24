clc; clear;
format rat
J1 = [5 0 3; 0 2 0; 3 0 5];

syms L real
factor(det(L*eye(3) - J1));

%L = 2
vl2 = [1; 0; -1];
vl22 = [0; 1; 0];
vl2a = [-1; 0; 1];
vl22a = [0; -1; 0];
(2*eye(3) - J1);
(2*eye(3) - J1)*vl22a;

%L = 8
vl8 = [1; 0; 1];
vl88 = [-1; 0; -1];
(8*eye(3) - J1);
(8*eye(3) - J1)*vl88;

vl2u = vl2/norm(vl2);
vl8u = vl8/norm(vl8);

cross(vl2u,vl8u);
R12 = [cross(vl2u,vl8u) vl2u vl8u];

[V,D] = eig(J1);

a = 2;
t = [0:0.01:4];
wx = a*sin(2*t);
wy = a*cos(2*t) - a;

plot(t,wx)
hold on
plot(t,wy)
