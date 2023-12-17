function [delta_x_plus, P_plus, NEES, NIS] = LKF(delta_x0, P0, delta_y, lmks_in_view_sim, F, Q, OMEGA, H_all, R, deltaX_true)
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
    H{i} = H_all{i}(idxs_H,:);
end

% preallocate vectors
delta_x_plus = zeros(n,npts_obs);
P_plus = zeros(n,n,npts_obs);

% ititialize state and uncertainty
delta_x_plus(:,1) = delta_x0;
P_plus(:,:,1) = P0;
k = 1;
delta_x_minus = delta_x0;
P_minus = P0;
NEES = zeros(1,npts_obs-1);
NIS = zeros(1,npts_obs-1);
for i = 2:npts_obs
    for j = 1:10
        [delta_x_minus, P_minus] = LKF_dynamicPrediction(delta_x_minus, P_minus, F(:,:,k), Q, OMEGA);
        k = k+1;
%         invPyykp1 = inv(P_minus);
    end
    dy = zeros(2*sum(lmks_in_view_sim(:,i)),1);
    dy(1:2:100) = lmks_in_view_sim(:,i);
    dy(2:2:100) = lmks_in_view_sim(:,i);
    dy_vis = delta_y(dy==1,i);
    if length(dy_vis) ~= 0
        [delta_x_plus(:,i), P_plus(:,:,i), innov_plus, S_k] = LKF_measurementUpdate(delta_x_minus, P_minus, dy_vis, H{i}, R);
        delta_x_minus = delta_x_plus(:,i);
        P_minus = P_plus(:,:,i);
        invPkp1 = inv(P_plus(:,:,i));
        
        if(deltaX_true)
            NEES(i-1) = (deltaX_true(:,i) - delta_x_plus(:,i))'*invPkp1*(deltaX_true(:,i) - delta_x_plus(:,i));
            NIS(i-1) = innov_plus'*S_k*innov_plus; 
        end
    end

end


end

