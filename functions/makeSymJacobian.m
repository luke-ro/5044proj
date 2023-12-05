function [jac_sym] = makeSymJacobian()
%MAKESYMJACOBIAN Summary of this function goes here
%   Detailed explanation goes here
syms x1 x2 x3 x4 x5 x6
r3 = (sqrt(x1^2+x2^2+x3^2))^3;
X = [x1,x2,x3,x4,x5,x6];
Ftil = [x4,x5,x6,-mu/r3*x1,-mu/r3*x2,-mu/r3*x3];
temp = matlabFunction(jacobian(Ftil,X));
jac_sym = @(st)temp(st(1),st(2),st(3));
end

