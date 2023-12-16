function [] = dataLKF(t_obs, uv_stacked_data, lmks_in_view, uv_stacked_nom, Q_LKF, P0_LKF, F, H, OMEGA, X_nom_N, const)

R = diag([const.sig_uv^2, const.sig_uv^2]);

delta_y = uv_stacked_data - uv_stacked_nom;
delta_x0 = zeros(6,1);
[delta_x_LKF, P_LKF, ~,~] = LKF(delta_x0, P0_LKF, delta_y, lmks_in_view, F, Q_LKF, OMEGA, H, R, false);

figure
% plotFilterResults(0, t_obs, 0, delta_x_LKF, P_LKF, "LKF Perturbation State Estimate from Data Log", "$$\delta$$")
plotFilterResults(0, t_obs, 0, delta_x_LKF+X_nom_N, P_LKF, "LKF Estimate from Data Log", "")
end

