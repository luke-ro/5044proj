function [] = dataEKF(t_obs, uv_stacked_data, lmks_in_view, P0_LKF, Q_EKF, gamma, const, R_CtoN, pos_lmks_A, pos_lmks_N)

R = diag([const.sig_uv^2, const.sig_uv^2]);

X0_true = [const.r0_nom_N; const.v0_nom_N];
[xhat_k_plus_hist, P_k_plus_hist, ~,~] = EKF(t_obs, X0_true, P0_LKF, const, gamma, Q_EKF, R_CtoN, pos_lmks_A, pos_lmks_N, uv_stacked_data, lmks_in_view, R, false);

figure
plotFilterResults(0, t_obs, 0, xhat_k_plus_hist, P_k_plus_hist, "EKF Estimate from Data Log","")

end

