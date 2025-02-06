function dX = SprottB(t,X)

% equations in numbers
R1 = 32e3;
R2 = 210;
R3 = 83e3;
R4 = 7.07;
R5 = 1e3;
R6 = 120;
R7 = 10e3;
R8 = 1e3;
L = 360e-6;
RL = 2.8;

C1 = 10e-9;
C2 = 1e-9;
V = 1.3;

mu = R4*R7/R6;

a1 = mu/10/R1/C1;
a2 = 1/L;
a3 = -(R4 + RL)/L;
a4 = V/C2/R3;
a5 = - mu/10/R2/C2;

dX = X;
x = X(1); y = X(2); z = X(3);
dX(1) = a1*y*z;
dX(2) = a2*x + a3*y;
dX(3) = a4 + a5*x*y;
end