function [delta_x_plus, delta_y, P_plus] = LKF(delta_x0, delta_y)
%LKF Linearized Kalman Filter
%   Detailed explanation goes here

n = length(delta_x0);
m = length(delta_y);
% preallocate vectors
delta_x_minus = zeros(n,m);
delta_x_plus = zeros(n,m);
P_minus = zeros(n,n,m);
P_plus = zeros(n,n,m);

% ititialize state and uncertainty
delta_x_plus(:,1) = x0_plus;
P_plus(:,:,1) = P0_plus;
% Kalman filter
for k = 1:m-1
   delta_x_minus(:,k+1) =  F*delta_x_plus(:,k);
   P_minus(:,:,k+1) = F*P_plus(:,:,k)*F' + Q;
   K = P_minus(:,:,k+1) *H'* inv(H*P_minus(:,:,k+1)*H' + R);
   delta_x_plus(:,k+1) = delta_x_minus(:,k+1) + K*(delta_y(:,k+1) - H*delta_x_minus(:,k+1)); 
   P_plus(:,:,k+1) = (eye(n,n) - K*H)*P_minus(:,:,k+1);

   P_minus(:,:,k+1) = 1/2*(P_minus(:,:,k+1) + P_minus(:,:,k+1)');
   P_plus(:,:,k+1) = 1/2*(P_plus(:,:,k+1) + P_plus(:,:,k+1)');
end


end

