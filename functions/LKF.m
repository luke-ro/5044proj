function [delta_x, delta_y, P] = LKF(delta_x0, P0, delta_y, F, Q, H, R)
%LKF Linearized Kalman Filter
%   Detailed explanation goes here

n = length(delta_x0);
m = length(delta_y);
% preallocate vectors
delta_x = zeros(n,m);
P = zeros(n,n,m);

% ititialize state and uncertainty
delta_x(:,1) = x0_plus;
P(:,:,1) = P0;

for k = 1:m-1
    [delta_x(:,k+1), P(:,:,k+1)] = LKFstep(delta_x(:,k), P(:,:,k), delta_y, F(:,:,k), Q, H{k}, R{k});
end

end

