function [xhat_k_plus,P_k_plus, innov_plus, S_k] = EKF_measurementUpdate(t, xhat_kplus1_minus, P_kplus1_minus, R_CtoN, pos_lmks_A, pos_lmks_N, const, nom_lmks_visible, i, y_true, R, n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [u_k, v_k, ~] = simMeasurements(t, xhat_kplus1_minus, R_CtoN(:,:,i), pos_lmks_A, const);
    H_tilde_kplus1 = zeros(100,6);
    lmk_idxs = 1:50;
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

    innov_plus = e_tilde_ykplus1;
    S_k = H_tilde_kplus1*P_kplus1_plus*H_tilde_kplus1' + R_mat;

end