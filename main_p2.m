% main_p2.m: part two main
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
close all
warning off

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
X0_nom_N = [const.r0_nom_N; const.v0_nom_N];

[X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(w_tilde, X0_nom_N, const);

% plot inertial orbit
title1 = "Simulated Noisy Trajectory in Inertial Frame";
plotOrbit(X_sim_N(1:3,:), title1)
figure
title2 = "Simulated Noisy States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,X_sim_N,title2,"")

% simulate measurements
[us, vs, sim_lmks_visible] = simMeasurements(t_obs, X_simObs_N, R_CtoN, pos_lmks_A, const);
uv_stacked_sim = stackUsVs(us,vs);
uv_stacked_nom = stackUsVs(u_nom,v_nom);

y_delta_sim = uv_stacked_sim - uv_stacked_nom;

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

% R = diag([const.sig_uv^2, const.sig_uv^2]);
% delta_X0 = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
delta_X0 = zeros(6,1);
% P0 = zeros(3,3);
% NOT properly accounting for time steps and visible landmarks together yet
% delta_y = (y_table_sim(:,3:4) - y_table_nom(:,3:4))';
% [delta_x_plus, P_plus] = LKF(delta_X0, P0, delta_y, bigF, Q, bigC, R);

% delta_x0 = X_deltaObs_N(:,1);
P0 = 0.5*eye(6);
R = diag([const.sig_uv^2, const.sig_uv^2]);
% Q = const.sig_w^2 * [const.Dt_int^3/3 0 0 const.Dt_int^2/2 0 0;
%                      0 const.Dt_int^3/3 0 0 const.Dt_int^2/2 0;
%                      0 0 const.Dt_int^3/3 0 0 const.Dt_int^2/2;
%                      const.Dt_int^2/2 0 0 const.Dt_int 0 0;
%                      0 const.Dt_int^2/2 0 0 const.Dt_int 0;
%                      0 0 const.Dt_int^2/2 0 0 const.Dt_int];

% Q = const.sig_w^2 * [600^3/3 0 0 600^2/2 0 0;
%                      0 600^3/3 0 0 600^2/2 0;
%                      0 0 600^3/3 0 0 600^2/2;
%                      600^2/2 0 0 600 0 0;
%                      0 600^2/2 0 0 600 0;
%                      0 0 600^2/2 0 0 600];

[delta_x_plus, P_plus, NEES, NIS] = LKF(delta_X0, P0, y_delta_sim, sim_lmks_visible, bigF, Q, bigC, R, X_nomObs_N);

figure
plotStates(t_obs,delta_x_plus,"","")

figure
hold on
plot(t_obs,delta_x_plus(1,:))
plot(t_obs,2*sqrt(squeeze(P_plus(1,1,:))))
plot(t_obs,-2*sqrt(squeeze(P_plus(1,1,:))))
hold off




Nsimruns = 50;
% NEESsamps = zeros(Nsimruns,length(t));
% NISsamps = zeros(Nsimruns,length(t));
% m0 = X_delta_N(:,1);
% P0 = eye(6);
% bigG = zeros(6,1,433);
% for ss=1:Nsimruns
%     %%%1. %%Generate true trajectory and measurements from system
%     
% %     xk_truehist = zeros(2,length(t));
%     ykhist = zeros(2,50,length(tvec));
%     xk = mvnrnd(m0,P0); %sample initial robot state
%     for jj=1:50 % go through all landmarks
%         
% %         wk = mvnrnd(zeros(1,2),Qtrue);
% %         xkp1 = F*xk' + G*u(jj) + wk';
% %         vkp1 = mvnrnd(zeros(1,1),Rtrue)';
% %         ykp1 = H*xkp1 + vkp1;
% %         
% %         xk_truehist(:,jj) = xkp1;
% %         ykhist(:,jj) = ykp1;
% %         xk = xkp1';
%         ykhist(1,jj,:) = us(jj,:);
%         ykhist(2,jj,:) = vs(jj,:);
%     end
% 
%     
%     %%% 2. Kalman Filter equations with simple NCV model
%     
%     %%Run the Kalman filter updates
%     mk = m0;
%     Pk = P0;
%     mk_filt_hist = zeros(2,length(t));
%     Pk_filt_hist = zeros(2,2,length(t));
%     innovk_hist = zeros(1,length(t));
%     Pyyk_hist = zeros(1,1,length(t)); %%store measurement innovation covar
%     NEESsshist = zeros(1,length(t));
%     NISsshist = zeros(1,length(t));
%     Rtrue = zeros(4);
%     for jj=1:length(tvec)
%         H = bigH{jj};
%         F = bigF(:,:,jj);
%         G = bigG(:,1,jj);
%         num_landmarks = sum(nom_lmks_visible(:,jj));
%         lmk_idxs = find(nom_lmks_visible(:,jj));
%         y_pred = nan(2,50);
%         for lm = 1:50
%         
%             %%Perform prediction step
%     %         mkp1_minus = F*mk' + G*u(jj);
%             mkp1_minus = F*mk;
%             Pkp1_minus = F*Pk*F' + Q;
%             
%             %%Compute Kalman gain
%             Pyykp1 = H*Pkp1_minus*H' + Rtrue;
%             Pyykp1 = 0.5*(Pyykp1 + Pyykp1');
%             Kkp1 = Pkp1_minus*H'/(H*Pkp1_minus*H' + Rtrue);
%             %%Perform measurement update step
%             ykp1_report = ykhist(:,:,jj); %simulate the reporting of data from sensor
%             ykp1_pred = H*mkp1_minus; %predicted measurement
%             y_pred(1,:) = ykp1_pred(1:2:length(ykp1_pred));
%             y_pred(2,:) = ykp1_pred(2:2:length(ykp1_pred));
%             innov_kp1 = ykp1_report - y_pred; %compute meas innovation
%             mkp1_plus = mkp1_minus + Kkp1*innov_kp1; %compute update to state mean
%             Pkp1_plus = (eye(6) - Kkp1*H)*Pkp1_minus; %compute update to covar
%             
%             mk = mkp1_plus';
%             mk_filt_hist(:,jj) = mkp1_plus;
%             Pk = Pkp1_plus;
%             Pk_filt_hist(:,:,jj)= Pkp1_plus;
%             innovk_hist(:,jj) = innov_kp1;
%             
%             %%Compute and store NEES and NIS statistics:
%             invPkp1 = inv(Pkp1_plus);
%             invPyykp1 = inv(Pyykp1);
%             NEESsshist(jj) = (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
%             NISsshist(jj) = innov_kp1'*invPyykp1*innov_kp1;    
%         end
%     end
%     NEESsamps(ss,:) = NEESsshist;
%     NISsamps(ss,:) = NISsshist;
%     
%     %%sanity check:
%     %Plot state estimation errors versus time
%     f=figure(5),
%     
%     subplot(211)
%     plot(t,xk_truehist(1,:)-mk_filt_hist(1,:),'Color',[0.5 0.5 0.5],'LineWidth',1), hold on
%     plot(t,2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'m','LineWidth',1)
%     plot(t,-2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'m','LineWidth',1)
%     ylabel('$\xi$ (m)','Interpreter','latex')
%     grid on, title('State Estimation Errors')
%     subplot(212)
%     plot(t,xk_truehist(2,:)-mk_filt_hist(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1), hold on
%     plot(t,2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'m','LineWidth',1)
%     plot(t,-2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'m','LineWidth',1)
%     grid on
%     xlabel('Time (s)')
%     ylabel('$\dot \xi$ (m)','Interpreter','latex')
% end





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

