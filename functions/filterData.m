function [] = filterData(t_obs, uv_stacked_data, lmks_in_view, uv_stacked_nom, Q_LKF, Q_EKF, P0_LKF, P0_EKF)

% LKF --------------------------------------------------------------------
delta_y = uv_stacked_data - uv_stacked_nom;
delta_x0 = zeros(6,1);
[delta_x_LKF, P_LKF, ~,~] = LKF(delta_x0, P0_LKF, delta_y, lmks_in_view, F, Q, OMEGA, H_all, R, deltaX_true);
% ------------------------------------------------------------------------

% EKF ---------------------------------------------------------------------

% -------------------------------------------------------------------------

end

