function [xhat_k_plus_hist, P_k_plus_hist, NEES, NIS] = EKF(t_obs, xhat0_plus, Phat0_plus, const, gamma, Q, R_CtoN, pos_lmks_A, pos_lmks_N, y_true, nom_lmks_visible, R, X_true)
%EKF Extended Kalman Filter
%   xhat0_plus : initial total state
%   Phat0_plus : initial covariance 

%   delta_y: cell array of stacked us and vs from simulated trajectorye
%   lmks_in_view: which landmarks are in the measurenemts delta_y
%   F: linearized dyanmics
%   Q:
%   H_all: for all us and vs regardless of if theyre in view or not (need
%       to get rid of extra landmark elemts
%   R:

n = length(xhat0_plus);
npts_obs = length(t_obs);
% npts_int = 4321;
steps = 1; % how many dynamic prediction steps we do at a time (for the ZOH noise inputs)


%% STEP 1
% Initialize at k = 0 the initial total state and covariance
% xhat0_plus at k = 0
% Phat0_plus at k = 0

xhat_k_plus = xhat0_plus;
P_k_plus = Phat0_plus;
% ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

% preallocate vectors
xhat_k_plus_hist = zeros(n,npts_obs);
P_k_plus_hist = zeros(n,n,npts_obs);
NEES = zeros(1,length(t_obs)-1);
NIS = zeros(1,length(t_obs)-1);

xhat_k_plus_hist(:,1) = xhat_k_plus;
P_k_plus_hist(:,:,1) = P_k_plus;
%% STEP 2

% lmk_idxs = 1:50;
for i = 2:npts_obs
    for j = 1:10
        [xhat_k_plus, P_k_plus] = EKF_dynamicPrediction(xhat_k_plus, P_k_plus, Q, gamma, const, steps);
    end
    xhat_kplus1_minus = xhat_k_plus;
    P_kplus1_minus = P_k_plus;

    [xhat_k_plus, P_k_plus, innov_plus, S_k] = EKF_measurementUpdate(t_obs(i), xhat_kplus1_minus, P_kplus1_minus, R_CtoN, pos_lmks_A, pos_lmks_N, const, nom_lmks_visible, i, y_true, R, n);
    xhat_k_plus_hist(:,i) = xhat_k_plus;
    P_k_plus_hist(:,:,i) = P_k_plus;

    invPkp1 = inv(P_k_plus);
    
    if(length(X_true)>100)
        NEES(i-1) = (X_true(:,i) - xhat_k_plus)'*invPkp1*(X_true(:,i) - xhat_k_plus);
        NIS(i-1) = innov_plus'*S_k*innov_plus; 
    end

end

% %calculate H from sim data
% H = cell(npts_obs,1);
% for i=1:length(H_all)
%     idxs_lmks = find(lmks_in_view_sim(:,i));
%     idxs_us = 2*(idxs_lmks-1)+1;
%     idxs_vs = 2*(idxs_lmks-1)+2;
%     idxs_H = sort([idxs_us;idxs_vs]);
%     H{i} = H_all{i}(idxs_H,:);
% end
% 
% % preallocate vectors
% delta_x_plus = zeros(n,npts_obs);
% P_plus = zeros(n,n,npts_obs);
% 
% % ititialize state and uncertainty
% delta_x_plus(:,1) = delta_x0;
% P_plus(:,:,1) = P0;
% k = 1;
% delta_x_minus = delta_x0;
% P_minus = P0;
% for i = 1:npts_obs-1
%     for j = 1:10
%         [delta_x_minus, P_minus] = LKF_dynamicPrediction(delta_x_minus, P_minus, F(:,:,k), Q(:,:,k));
%         k = k+1;
% %         invPyykp1 = inv(P_minus);
%     end
%     dy = zeros(2*sum(lmks_in_view_sim(:,i)),1);
%     dy(1:2:100) = lmks_in_view_sim(:,i);
%     dy(2:2:100) = lmks_in_view_sim(:,i);
%     dy_vis = delta_y(find(dy==1),i);
%     [delta_x_plus(:,i), P_plus(:,:,i), innov_plus, S_k] = LKF_measurementUpdate(delta_x_minus, P_minus, dy_vis, H{i}, R);
%     
%     invPkp1 = inv(P_plus(:,:,i));
%     NEES(i) = (x_nomObs_N(:,i) - delta_x_plus(:,i))'*invPkp1*(x_nomObs_N(:,i) - delta_x_plus(:,i));
%     P_minus = P_plus(:,:,i);
%     NIS(i) = innov_plus'*S_k*innov_plus; 
% end


end

