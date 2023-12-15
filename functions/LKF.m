function [delta_x_plus, P_plus] = LKF(delta_x0, P0, delta_y, lmks_in_view_sim, F, Q, H_all, R)
%LKF Linearized Kalman Filter
%   delta_y: cell array of stacked us and vs from simulated trajectorye
%   lmks_in_view: which landmarks are in the measurenemts delta_y
%   F: linearized dyanmics
%   Q:
%   H_all: for all us and vs regardless of if theyre in view or not (need
%       to get rid of extra landmark elemts
%   R:



n = length(delta_x0);
npts_obs = length(delta_y);
% npts_int = length(F);

%calculate H from sim data
H = cell(npts_obs,1);
for i=1:length(H_all)
    idxs_lmks = find(lmks_in_view_sim(:,i));
    idxs_us = 2*(idxs_lmks-1)+1;
    idxs_vs = 2*(idxs_lmks-1)+2;
    idxs_H = sort([idxs_us;idxs_vs]);
    i
    H{i} = H_all{i}(idxs_H,:);
end

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
    [delta_x_plus(:,i), P_plus(:,:,i)] = LKF_measurementUpdate(delta_x_minus, P_minus, delta_y(:,i), H(:,:,i), R);
end

% count = 1;
% for k = 1:npts_int-1
%     [delta_x_minus(:,k+1), P_minus(:,:,k+1)] = LKF_dynamicPrediction(delta_x(:,k), P_minus(:,:,k), F(:,:,k), Q(:,:,k));
%     if mod(npts_int/npts_obs, k) == 0
%         [delta_x_plus(:,count), P_plus(:,:,count)] = LKF_measurementUpdate(delta_x_minus(:,k), P_minus(:,:,k), delta_y(:,count), H(:,:,count), R);
%         count = count + 1;
%     end
% end

end

