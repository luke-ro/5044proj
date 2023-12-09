function [delta_x_plus, P_plus] = LKFstep(delta_x_plus0, P_plus0, delta_y, F, Q, H, R)
%LKFSTEP One timestep of the linearized kalman filter
%   Detailed explanation goes here

n = length(delta_x0);
m = length(delta_y);

for k = 1:m-1
   delta_x_minus =  F*delta_x_plus0;
   P_minus = F*P_plus0*F' + Q;
   K = P_minus *H'* inv(H*P_minus*H' + R);
   delta_x_plus = delta_x_minus + K*(delta_y - H*delta_x_minus); 
   P_plus = (eye(n,n) - K*H)*P_minus;

%    P_minus = 1/2*(P_minus + P_minus');
   P_plus = 1/2*(P_plus + P_plus');
end

end

