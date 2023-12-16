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
% C_w_tilde = diag([zeros(1,3), (const.sig_w^2)*ones(1,3)]);
C_w_tilde = diag([zeros(1,3), (const.sig_w^2)*ones(1,3)]);
% C_w_tilde = diag([(const.sig_w^2)*ones(1,3), zeros(1,3)]);
w_tilde = mvnrnd(zeros(1,n), C_w_tilde, npts_int)';
[X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(w_tilde, X0_nom_N, const);

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

% W = diag([const.sig_w^2,const.sig_w^2,const.sig_w^2]);
% Q = zeros(size(bigA));
% for i = 1:length(bigA)
%     eZ = expm(const.Dt_int*[-bigA(:,:,i), gamma*W*gamma'; zeros(n,n), bigA(:,:,i)']);
%     Q(:,:,i) = eZ(n+1:end, n+1:end)' * eZ(1:n, n+1:end);
% %     w = mvnrnd(zeros(1,n), Q, npts_int);  
% end

delta_X0 = deltaX_sim_N(:,1);
P0 = blkdiag(0.01*eye(3), 1e-6*eye(3));
R = diag([const.sig_uv^2, const.sig_uv^2]);
Q = const.sig_w^2 * [const.Dt_int^3/3 0 0 const.Dt_int^2/2 0 0;
                     0 const.Dt_int^3/3 0 0 const.Dt_int^2/2 0;
                     0 0 const.Dt_int^3/3 0 0 const.Dt_int^2/2;
                     const.Dt_int^2/2 0 0 const.Dt_int 0 0;
                     0 const.Dt_int^2/2 0 0 const.Dt_int 0;
                     0 0 const.Dt_int^2/2 0 0 const.Dt_int];

% [delta_x_plus, P_plus, NEES, NIS] = LKF(delta_X0, P0, y_delta_sim, sim_lmks_visible, bigF, Q, bigC, R, X_nomObs_N);
[delta_x_plus, P_plus, NEES, NIS] = LKF(delta_X0, P0, y_delta_sim, nom_lmks_visible, bigF, Q, OMEGA, bigC, R, X_nomObs_N);

figure
plotStates(t_obs,delta_x_plus,"","")

figure
hold on
plot(t, deltaX_sim_N(1,:))
plot(t_obs,delta_x_plus(1,:))
plot(t_obs,deltaX_sim_N(1,1:10:end)+2*sqrt(squeeze(P_plus(1,1,:))'))
plot(t_obs,deltaX_sim_N(1,1:10:end)-2*sqrt(squeeze(P_plus(1,1,:))'))
legend("True $\delta x$","LKF $\delta x$", "$+2\sigma$","$+2\sigma$","Interpreter", "Latex")
% ylim([-20,20])
hold off

figure
hold on
plot(t, deltaX_sim_N(2,:))
plot(t_obs,delta_x_plus(2,:))
plot(t_obs,deltaX_sim_N(2,1:10:end)+2*sqrt(squeeze(P_plus(2,2,:))'))
plot(t_obs,deltaX_sim_N(2,1:10:end)-2*sqrt(squeeze(P_plus(2,2,:))'))
legend("True $\delta x$","LKF $\delta x$", "$+2\sigma$","$+2\sigma$","Interpreter", "Latex")
% ylim([-20,20])
hold off

figure
hold on
plot(t, deltaX_sim_N(3,:))
plot(t_obs,delta_x_plus(3,:))
plot(t_obs,deltaX_sim_N(3,1:10:end)+2*sqrt(squeeze(P_plus(3,3,:))'))
plot(t_obs,deltaX_sim_N(3,1:10:end)-2*sqrt(squeeze(P_plus(3,3,:))'))
legend("True $\delta x$","LKF $\delta x$", "$+2\sigma$","$+2\sigma$","Interpreter", "Latex")
% ylim([-20,20])
hold off




Nsimruns = 50;




%%DO NEES TEST:
% epsNEESbar = mean(NEESsamps,1);
epsNEESbar = mean(NEES,1);
alphaNEES = 0.05; %%significance level
Nnx = Nsimruns*length(bigF(:,:,1));
%%compute intervals:
r1x = chi2inv(alphaNEES/2, Nnx )./ Nsimruns
r2x = chi2inv(1-alphaNEES/2, Nnx )./ Nsimruns

figure()
% title("NEES Estimation")
title('NEES Estimation Results','FontSize',14)
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon}_x$','Interpreter','latex', 'FontSize',14)
xlabel('time step, k','FontSize',14)
% title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound'),grid on
% saveas(f,'NEESResults.png','png')
hold off

%%DO NIS TEST:
epsNISbar = mean(NIS,1);
alphaNIS = 0.05; 

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
% title("NIS Estimation")
for i=1:length(bigC)
    Nny = Nsimruns*size(bigC{i},1);
    %%compute intervals:
    r1y = chi2inv(alphaNIS/2, Nny )./ Nsimruns;
    r2y = chi2inv(1-alphaNIS/2, Nny )./ Nsimruns;
    
    
    plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
    plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
end
ylabel('NIS statistic, $\bar{\epsilon}_y$','Interpreter','latex','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound'),grid on
% saveas(f,'NISResults.png','png')

