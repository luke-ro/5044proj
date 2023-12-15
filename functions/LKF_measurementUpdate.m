function [delta_x_plus, P_plus, innov_plus, S_k] = LKF_measurementUpdate(delta_x_minus, P_minus, delta_y, H, R)
%LKFSTEP One timestep of the linearized kalman filter
%   Detailed explanation goes here

n = length(delta_x_minus);
m = length(delta_y);

R_single = R;
for i = 1:size(H,1)/2-1
    R = blkdiag(R, R_single);
end
% invTerm = (H*P_minus*H' + R)\eye(size(H*P_minus*H' + R));
% invTerm = inv(H*P_minus*H' + R);
% K = (P_minus * H') * invTerm;
K = (P_minus * H') * inv(H*P_minus*H' + R);
innov_plus = delta_y - H*delta_x_minus;
% innov_plus = delta_y;
delta_x_plus = delta_x_minus + K*innov_plus; 
P_plus = (eye(n,n) - K*H)*P_minus;

% ensure it stays positive definite
P_plus = 1/2*(P_plus + P_plus');

S_k = H*P_minus*H' + R;
end

