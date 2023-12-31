function [J] = dyn_jacobian(x,const)
x1 = x(1);
x2 = x(2);
x3 = x(3);
% x4 = X(4);
% x5 = X(5);
% x6 = X(6);

Mua = const.mu;
denom = (x1^2+x2^2+x3^2)^(5/2);
denom1 = (x1^2+x2^2+x3^2)^(3/2);
J = zeros(6,6);
% J(4,1) = -Mua*((2*x1^2-x2^2-x3^2)/denom);
J(4,1) = -Mua*((1/denom1)-(3*x1^2/denom));
J(5,1) = 3*Mua*x1*x2/denom;
J(6,1) = 3*Mua*x1*x3/denom;

J(4,2) = 3*Mua*x1*x2/denom;
% J(5,2) = -Mua*((2*x2^2-x1^2-x3^2)/denom);
J(5,2) = -Mua*((1/denom1)-(3*x2^2/denom));
J(6,2) = 3*Mua*x2*x3/denom;

J(4,3) = 3*Mua*x1*x3/denom;
J(5,3) = 3*Mua*x2*x3/denom;
% J(6,3) = -Mua*((2*x3^2-x1^2-x2^2)/denom);
J(6,3) = -Mua*((1/denom1)-(3*x3^2/denom));

J(1:3,4:6) = eye(3);
end

