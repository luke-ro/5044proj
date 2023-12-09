function [func] = makeSymH(constants)
%MAKESYMH Summary of this function goes here
%   Detailed explanation goes here
l = sym("l",[3 1]);
% r = sym("r",[3,1]);
% rdot = sym("rdot",[3 1]);
% x = [r;rdot];
x = sym("x",[6,1]);
ic = sym("ic",[3,1]);
jc = sym("jc",[3,1]);
kc = sym("kc",[3,1]);
uv = [(constants.f*dot((l-x(1:3)),ic)/(dot((l-x(1:3)),kc))) + constants.u0;
      (constants.f*dot((l-x(1:3)),jc)/(dot((l-x(1:3)),kc))) + constants.v0];

symH = matlabFunction(jacobian(uv,x));

func = @(x0,lmk,R) symH(R(1,1),R(2,1),R(3,1),R(1,2),R(2,2),R(3,2),R(1,3),R(2,3),R(3,3),lmk(1),lmk(2),lmk(3),x0(1),x0(2),x0(3));
end

