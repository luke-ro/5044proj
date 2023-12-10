function [delta_x_plus, P_plus] = LKF(delta_x0, P0, delta_y, F, Q, H, R)
%LKF Linearized Kalman Filter
%   Detailed explanation goes here

n = length(delta_x0);
npts_obs = length(delta_y);
% npts_int = length(F);

% preallocate vectors
delta_x_plus = zeros(n,npts_obs);
P_plus = zeros(n,n,npts_obs);

% ititialize state and uncertainty
delta_x_plus(:,1) = delta_x0;
P_plus(:,:,1) = P0;
k = 1;
for i = 1:npts_obs-1
    for j = 1:10
        [delta_x_minus, P_minus] = LKF_dynamicPrediction(delta_x_minus, P_minus, F(:,:,k), Q(:,:,k));
        k = k+1;
    end
    [delta_x_plus(:,i), P_plus(:,:,i)] = LKF_measurementUpdate(delta_x_minus, P_minus, delta_y(:,i), H(:,:,i), R(:,:,i));
end

% count = 1;
% for k = 1:npts_int-1
%     [delta_x_minus(:,k+1), P_minus(:,:,k+1)] = LKF_dynamicPrediction(delta_x(:,k), P_minus(:,:,k), F(:,:,k), Q(:,:,k));
%     if mod(npts_int/npts_obs, k) == 0
%         [delta_x_plus(:,count), P_plus(:,:,count)] = LKF_measurementUpdate(delta_x_minus(:,k), P_minus(:,:,k), delta_y(:,count), H(:,:,count), R(:,:,count));
%         count = count + 1;
%     end
% end

end

