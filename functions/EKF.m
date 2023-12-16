function [xhat_k_plus_hist, P_k_plus_hist] = EKF(xhat0_plus, Phat0_plus, const, gamma, Q, R_CtoN, pos_lmks_A, pos_lmks_N, y_true, nom_lmks_visible, R)
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
npts_obs = 432;
npts_int = 4321;


%% STEP 1
% Initialize at k = 0 the initial total state and covariance
% xhat0_plus at k = 0
% Phat0_plus at k = 0

xhat_k_plus = xhat0_plus;
P_k_plus = Phat0_plus;
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

% preallocate vectors
xhat_k_plus_hist = zeros(n,npts_obs+1);
P_k_plus_hist = zeros(n,n,npts_obs+1);

xhat_k_plus_hist(:,1) = xhat_k_plus;
P_k_plus_hist(:,:,1) = P_k_plus;
%% STEP 2

k = 1;
lmk_idxs = 1:50;
for i = 1:npts_obs - 1
    for j = 1:10
        ode_fun = @(t,X) dynamics(t,X,const,zeros(n,npts_int));
        [t,X_sim] = ode45(ode_fun,[0 60],xhat_k_plus,ode_options);
    %     xhat_kplus1_minus = X_sim_N(:,1);
        xhat_k_plus = X_sim(end,:)';
    
        %calcuate jacobian accoring to nominal traj
        Atilde = dyn_jacobian(xhat_k_plus,const);
        
        %calcuate f using matrix expm
        Ahat = [Atilde zeros(6,1); zeros(1,7)];
        phim_hat = expm(Ahat*const.Dt_int);
        Ftilde = phim_hat(1:6, 1:6);
        Gtilde = phim_hat(1:6, 7);
        omega = const.Dt_int * gamma;
    
        P_kplus1_minus = Ftilde*P_k_plus*Ftilde' + omega*Q*omega';
        P_k_plus = P_kplus1_minus;
        k = k+1;
    end
    xhat_kplus1_minus = xhat_k_plus;

    [u_k, v_k, ~] = simMeasurements(const.Dt_int, xhat_kplus1_minus, R_CtoN(:,:,i), pos_lmks_A, const);
    H_tilde_kplus1 = zeros(100,6);
    for ii = 1:2:100
        H_tilde_kplus1(ii:ii+1,:) = dyn_jacobian_H(xhat_kplus1_minus, pos_lmks_N(:,lmk_idxs((ii+1)/2),i), R_CtoN(:,:,i), const);
    end
    u_k = u_k(nom_lmks_visible(:,i));
    v_k = v_k(nom_lmks_visible(:,i));
    yhat_kplus1 = zeros(2*sum(nom_lmks_visible(:,i)),1);
    yhat_kplus1(1:2:length(u_k)*2) = u_k;
    yhat_kplus1(2:2:length(v_k)*2) = v_k;

    true_mask = zeros(2*sum(nom_lmks_visible(:,i)),1);
    true_mask(1:2:100) = nom_lmks_visible(:,i);
    true_mask(2:2:100) = nom_lmks_visible(:,i);
    true_mask = logical(true_mask);
    e_tilde_ykplus1 = y_true(true_mask,i) - yhat_kplus1;

    H_tilde_kplus1 = H_tilde_kplus1(true_mask,:);
    R_single = R;
    R_mat = R;
    for jj = 1:size(H_tilde_kplus1,1)/2-1
        R_mat = blkdiag(R_mat, R_single);
    end

    invTerm = (H_tilde_kplus1 * P_kplus1_minus * H_tilde_kplus1' + R_mat)\eye(size((H_tilde_kplus1 * P_kplus1_minus * H_tilde_kplus1' + R_mat)));
    K_tilde_kplus1 = P_kplus1_minus * H_tilde_kplus1' * invTerm;

    xhat_kplus1_plus = xhat_kplus1_minus + K_tilde_kplus1*e_tilde_ykplus1;
    P_kplus1_plus = (eye(n) - K_tilde_kplus1*H_tilde_kplus1)*P_kplus1_minus;
    
    xhat_k_plus = xhat_kplus1_plus;
    P_k_plus = P_kplus1_plus;

    xhat_k_plus_hist(:,i+1) = xhat_k_plus;
    P_k_plus_hist(:,:,i+1) = P_k_plus;
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

