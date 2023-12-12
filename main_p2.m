% main_p2.m: part two main
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
close all

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
% C_w_tilde = diag([zeros(1,3), (const.sig_w^2)*ones(1,3)]);
C_w_tilde = diag([(const.sig_w^2)*ones(1,3), zeros(1,3)]);
w_tilde = mvnrnd(zeros(1,n), C_w_tilde, npts_int)';

[X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(w_tilde, const);

% plot inertial orbit
title = "Simulated Noisy Trajectory in Inertial Frame";
plotOrbit(X_sim_N(1:3,:), title)
figure
title = "Simulated Noisy States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,X_sim_N,title,"")

% simulate measurements
[us, vs, sim_lmks_visible] = simMeasurements(t_obs, X_simObs_N, R_CtoN, pos_lmks_A, const);

% generate a y table with the simulated data
y_table_sim = genYTable(t_obs,us,vs,(sim_lmks_visible & nom_lmks_visible));
y_table_nom = genYTable(t_obs,u_nom,v_nom,(sim_lmks_visible & nom_lmks_visible));

plotMeasurements(t_obs, us, vs, sim_lmks_visible, 1:10, "Simulated noisy measurements")

%% DT simulation

gamma = zeros(6,3);
gamma(4:6,:) = eye(3);
W = diag([const.sig_w^2,const.sig_w^2,const.sig_w^2]);
Q = zeros(size(bigA));
for i = 1:length(bigA)
    eZ = expm(const.Dt_int*[-bigA(:,:,i), gamma*W*gamma'; zeros(n,n), bigA(:,:,i)']);
    Q(:,:,i) = eZ(n+1:end, n+1:end)' * eZ(1:n, n+1:end);
%     w = mvnrnd(zeros(1,n), Q, npts_int);  
end

R = diag([const.sig_uv^2, const.sig_uv^2]);
delta_X0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
P0 = zeros(3,3);
% NOT properly accounting for time steps and visible landmarks together yet
% delta_y = (y_table_sim(:,3:4) - y_table_nom(:,3:4))';
% [delta_x_plus, P_plus] = LKF(delta_X0, P0, delta_y, bigF, Q, bigC, R);



