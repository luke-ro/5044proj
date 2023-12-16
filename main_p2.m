% main_p2.m: part two main
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
close all
% warning off
rng(100)

addpath('./functions/');
addpath('./Data/');
load('orbitdetermination-finalproj_data_2023_11_14.mat')
load('P1_vars.mat')

n = 6;
m = 2;
% npts_int = length(0:const.Dt_int:const.tf_int);
npts_int = length(0:const.Dt_int:const.tf_obs);
npts_obs = length(0:const.Dt_obs:const.tf_obs);

%% Simulate full nonlinear dynamics with process noise
X0_nom_N = [const.r0_nom_N; const.v0_nom_N];

% simulate with process noise ---------------------------------------------
X0_delta = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
C_w_tilde = diag([zeros(1,3), (const.sig_w^2)*ones(1,3)]);
w_tilde = mvnrnd(zeros(1,n), C_w_tilde, npts_int)';
[X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(w_tilde, X0_nom_N+X0_delta, const);

% simulate measurements
[us, vs, sim_lmks_visible] = simMeasurements(t_obs, X_simObs_N, R_CtoN, pos_lmks_A, const);
uv_stacked_sim = stackUsVs(us,vs);
uv_stacked_nom = stackUsVs(u_nom,v_nom);
% -------------------------------------------------------------------------

% % test with no process noise ----------------------------------------------
% X0_delta = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
% X0 = X0_nom_N + X0_delta;
% [X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(zeros(6,6,npts_int), X0, const);
% % simulate measurements
% [us, vs, sim_lmks_visible] = simMeasurements(t_obs, X_simObs_N, R_CtoN, pos_lmks_A, const);
% uv_stacked_sim = stackUsVs(us,vs);
% uv_stacked_nom = stackUsVs(u_nom,v_nom);
% % -------------------------------------------------------------------------

deltaX_sim_N = X_sim_N - X_nom_N;
y_delta_sim = uv_stacked_sim - uv_stacked_nom;

% generate a y table with the simulated data
y_table_sim = genYTable(t_obs,us,vs,(sim_lmks_visible & nom_lmks_visible));
y_table_nom = genYTable(t_obs,u_nom,v_nom,(sim_lmks_visible & nom_lmks_visible));

% plot inertial orbit
title1 = "Simulated Noisy Trajectory in Inertial Frame";
plotOrbit(X_sim_N(1:3,:), title1)
figure
% plot states
title2 = "Simulated Noisy States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,X_sim_N,title2,"")
figure
title2 = "Simulated Noisy Perturbation States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,deltaX_sim_N,title2,"$$\delta$$")

plotMeasurements(t_obs, us, vs, sim_lmks_visible, 1:10, "Simulated noisy measurements")

%% DT simulation

gamma = zeros(6,6);
gamma(4:6,4:6) = eye(3);
OMEGA = const.Dt_int*gamma;
R = diag([const.sig_uv^2, const.sig_uv^2]);

delta_X0 = deltaX_sim_N(:,1);
P0_LKF = blkdiag(0.01*eye(3), 1e-6*eye(3));
Q_LKF = 1000*const.sig_w^2 * [const.Dt_int^3/3 0 0 const.Dt_int^2/2 0 0;
                     0 const.Dt_int^3/3 0 0 const.Dt_int^2/2 0;
                     0 0 const.Dt_int^3/3 0 0 const.Dt_int^2/2;
                     const.Dt_int^2/2 0 0 const.Dt_int 0 0;
                     0 const.Dt_int^2/2 0 0 const.Dt_int 0;
                     0 0 const.Dt_int^2/2 0 0 const.Dt_int];


[delta_x_plus, P_plus, NEES, NIS] = LKF(delta_X0, P0_LKF, y_delta_sim, nom_lmks_visible, bigF, Q_LKF, OMEGA, bigC, R, deltaX_sim_N);
figure
plotFilterResults(t, t_obs, deltaX_sim_N, delta_x_plus, P_plus, "LKF Typical Run", "$$\delta$$")

P0_EKF = blkdiag(0.01*eye(3), 1e-6*eye(3));
Q_EKF = 1000*const.sig_w^2 * [const.Dt_int^3/3 0 0 const.Dt_int^2/2 0 0;
                     0 const.Dt_int^3/3 0 0 const.Dt_int^2/2 0;
                     0 0 const.Dt_int^3/3 0 0 const.Dt_int^2/2;
                     const.Dt_int^2/2 0 0 const.Dt_int 0 0;
                     0 const.Dt_int^2/2 0 0 const.Dt_int 0;
                     0 0 const.Dt_int^2/2 0 0 const.Dt_int];
X0_true = [const.r0_nom_N; const.v0_nom_N];
[xhat_k_plus_hist, P_k_plus_hist, NEES_EKF, NIS_EKF] = EKF(t_obs, X0_true, P0_EKF, const, gamma, Q_EKF, R_CtoN, pos_lmks_A, pos_lmks_N, uv_stacked_sim, nom_lmks_visible, R, X_simObs_N);
figure
plotFilterResults(t, t_obs, X_sim_N, xhat_k_plus_hist, P_k_plus_hist, "EKF Results","")


Nsimruns = 10;
% X0_delta = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
P0 = blkdiag((1e-5)^2*eye(3), (1e-7)^2*eye(3));
calcNEESNIS(Nsimruns,t_obs,P0_LKF,P0,C_w_tilde,Q_LKF,R,OMEGA,0.05,0.05, const)
calcEKFNEESNIS(Nsimruns,t_obs,P0_EKF,P0,C_w_tilde, Q_EKF, R, 0.05, 0.05, const, gamma)


for i = 1:433
    if any(diag([P_k_plus_hist(:,:,i)]) < 0)
        fprintf("Negative diagonal elements at k = %d \n", i)
    end
end

[Y,lmks_vis] = parseYtable(y_table);
dataEKF(t_obs, Y, lmks_vis, P0_EKF, Q_EKF, gamma, const, R_CtoN, pos_lmks_A, pos_lmks_N)
dataLKF(t_obs, Y, lmks_vis, uv_stacked_nom, Q_LKF, P0_LKF, bigF, bigC, OMEGA, X_nomObs_N, const)

