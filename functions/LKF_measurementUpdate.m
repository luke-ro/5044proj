function [delta_x_plus, P_plus] = LKF_measurementUpdate(delta_x_minus, P_minus, delta_y, H, R)
%LKFSTEP One timestep of the linearized kalman filter
%   Detailed explanation goes here

n = length(delta_x0);
m = length(delta_y);

K = P_minus *H'* inv(H*P_minus*H' + R);
delta_x_plus = delta_x_minus + K*(delta_y - H*delta_x_minus); 
P_plus = (eye(n,n) - K*H)*P_minus;

% ensure it stays positive definite
P_plus = 1/2*(P_plus + P_plus');

end

