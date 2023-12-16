function [outputArg1,outputArg2] = calcNEESNIS(n_runs, t_obs, P0_LKF, P0, C_w_tilde, Qkf, R, OMEGA, alpha_nees, alpha_nis, const)
%CALCNEES Summary of this function goes here
%   Detailed explanation goes here

load("P1_vars.mat","bigF","bigC","u_nom","v_nom","nom_lmks_visible","X_nomObs_N")
load('data/orbitdetermination-finalproj_data_2023_11_14.mat',"R_CtoN","pos_lmks_A")

X0_nom = [const.r0_nom_N;const.v0_nom_N];
%calculate nominal traj and measurements
% [X_nom_N, t, X_nomObs_N, t_obs]= simNLdynamics(w_tilde, X0_nom, const);
% [us, vs, nom_lmks_visible] = simMeasurements(t_obs, X_nomObs_N, R_CtoN, pos_lmks_A, const);
Y_nom_N = stackUsVs(u_nom,v_nom);


NEES_hist = zeros(n_runs, length(t_obs)-1);
NIS_hist = zeros(n_runs,length(t_obs)-1);

for i = 1:n_runs
    % calculate initial state using m0 and P0:
    delta_X0_rand = mvnrnd(zeros(6,1),P0)';
    w_tilde = mvnrnd(zeros(1,6), C_w_tilde, length(0:const.Dt_int:const.tf_int))';

    %calculate true trajectory and measurements using non-linear dynamics
    [~, ~, X_simObs_N, t_obs] = simNLdynamics(w_tilde, X0_nom+delta_X0_rand, const);

    [us_sim, vs_sim, sim_lmks_visible] = simMeasurements(t_obs, X_simObs_N, R_CtoN, pos_lmks_A, const);
%     uv_stacked_sim = stackUsVs(us,vs);
    Y_sim_N = stackUsVs(us_sim,vs_sim);
    
    Y_delta = Y_sim_N - Y_nom_N;

    % run the kalman filter 
    % X_est: (n x time) matrix of estimated states
    % P_yy:  (p x p X time) matrix of Pyy=0.5*H_kp1*P_minus_kp1*H_kp1'
    % P_plus: (n x n x time) matrix of covariance matrix after meas update
    % y_innov: (p x n x time) matrix of dynamics predicted y - meas y
    deltaX_true = X_simObs_N-X_nomObs_N;
    
    [delta_x_plus, P_plus, NEES, NIS]= LKF(delta_X0_rand, P0_LKF, Y_delta, nom_lmks_visible, bigF, Qkf, OMEGA, bigC, R, deltaX_true);

    NEES_hist(i,:) = NEES;
    NIS_hist(i,:) = NIS;
end

epsNEESbar = mean(NEES_hist,1);
% alphaNEES = 0.05; %%significance level
Nnx = n_runs*length(bigF(:,:,1));
%%compute intervals:
r1x = chi2inv(alpha_nees/2, Nnx ) ./ n_runs;
r2x = chi2inv(1-alpha_nees/2, Nnx ) ./ n_runs;

figure()
% title("NEES Estimation")
% title('LKF NEES Estimation Results','FontSize',14)
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
xlabel('time step, k','FontSize',14)
title('LKF - NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on
% saveas(f,'NEESResults.png','png')

epsNISbar = mean(NIS,1);

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
title("LKF NIS Estimation")
for i=1:length(bigC)
    Nny = n_runs*size(bigC{i},1);
    %%compute intervals:
    r1y = chi2inv(alpha_nis/2, Nny )./ n_runs;
    r2y = chi2inv(1-alpha_nis/2, Nny )./ n_runs;
    
    
    plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
    plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
end
ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
xlabel('time step, k','FontSize',14)
title('LKF - NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on
% saveas(f,'NISResults.png','png')
end

