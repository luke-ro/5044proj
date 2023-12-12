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
% npts = length(0:const.Dt_int:const.tf_int);
npts = length(0:const.Dt_obs:const.tf_obs);

%% Simulate full nonlinear dynamics with process noise
C_w_tilde = diag([zeros(1,3), (const.sig_w^2)*ones(1,3)]);
w_tilde = mvnrnd(zeros(1,n), C_w_tilde, npts)';

[X_sim_N, t] = simNLdynamics(w_tilde, const);

% plot inertial orbit
title = "Simulated Noisy Trajectory in Inertial Frame";
plotOrbit(X_sim_N(1:3,:), title)
figure
title = "Simulated Noisy States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,X_sim_N,title,"")

tvec = t(1):600:t(end);
% simulate measurements
[us, vs, lmks_visible] = simMeasurements(tvec, X_sim_N, R_CtoN, pos_lmks_A, const);

% generate a y table with the simulated data
y_table_sim = genYTable(tvec,us,vs,lmks_visible);

plotMeasurements(tvec, us, vs, lmks_visible, 1:10, "Simulated noisy measurements")

%% DT simulation

gamma = zeros(6,3);
gamma(4:6,:) = eye(3);
W = diag([const.sig_w^2,const.sig_w^2,const.sig_w^2]);
for i = 1:length(bigA)
    eZ = expm(const.Dt_int*[-bigA(:,:,i), gamma*W*gamma'; zeros(n,n), bigA(:,:,i)']);
    Q = eZ(n+1:end, n+1:end)' * eZ(1:n, n+1:end);
    w = mvnrnd(zeros(1,n), Q, npts);  
end

[us, vs, lmks_visible] = simMeasurements(t, X_sim_N, R_CtoN, pos_lmks_A, const);
plotMeasurements(t, us, vs, lmks_visible, 10, "Full Nonlinear Measurement Simulation")

% R = diag([sig_uv^2, sig_uv^2]);
% v_u = mvnrnd(zeros(1,n), Q, npts);
% v_v = mvnrnd(zeros(1,n), R, npts);



Nsimruns = 50;
NEESsamps = zeros(Nsimruns,length(t));
NISsamps = zeros(Nsimruns,length(t));
m0 = X_delta_N(:,1);
P0 = eye(6);
bigG = zeros(6,1,433);
for ss=1:Nsimruns
    %%%1. %%Generate true trajectory and measurements from system
    
%     xk_truehist = zeros(2,length(t));
    ykhist = zeros(2,50,length(tvec));
    xk = mvnrnd(m0,P0); %sample initial robot state
    for jj=1:50 % go through all landmarks
        
%         wk = mvnrnd(zeros(1,2),Qtrue);
%         xkp1 = F*xk' + G*u(jj) + wk';
%         vkp1 = mvnrnd(zeros(1,1),Rtrue)';
%         ykp1 = H*xkp1 + vkp1;
%         
%         xk_truehist(:,jj) = xkp1;
%         ykhist(:,jj) = ykp1;
%         xk = xkp1';
        ykhist(1,jj,:) = us(jj,:);
        ykhist(2,jj,:) = vs(jj,:);
    end

    
    %%% 2. Kalman Filter equations with simple NCV model
    
    %%Run the Kalman filter updates
    mk = m0;
    Pk = P0;
    mk_filt_hist = zeros(2,length(t));
    Pk_filt_hist = zeros(2,2,length(t));
    innovk_hist = zeros(1,length(t));
    Pyyk_hist = zeros(1,1,length(t)); %%store measurement innovation covar
    NEESsshist = zeros(1,length(t));
    NISsshist = zeros(1,length(t));
    Rtrue = zeros(48);
    for jj=1:length(tvec)
        H = bigC{jj};
        F = bigF(:,:,jj);
        G = bigG(:,1,jj);
        for lm = 1:50
        
            %%Perform prediction step
    %         mkp1_minus = F*mk' + G*u(jj);
            mkp1_minus = F*mk;
            Pkp1_minus = F*Pk*F' + Q;
            
            %%Compute Kalman gain
            Pyykp1 = H*Pkp1_minus*H' + Rtrue;
            Pyykp1 = 0.5*(Pyykp1 + Pyykp1');
            Kkp1 = Pkp1_minus*H'/(H*Pkp1_minus*H' + Rtrue);
            %%Perform measurement update step
            ykp1_report = ykhist(:,:,jj); %simulate the reporting of data from sensor
            ykp1_pred = H*mkp1_minus; %predicted measurement
            innov_kp1 = ykp1_report - ykp1_pred; %compute meas innovation
            mkp1_plus = mkp1_minus + Kkp1*innov_kp1; %compute update to state mean
            Pkp1_plus = (eye(6) - Kkp1*H)*Pkp1_minus; %compute update to covar
            
            mk = mkp1_plus';
            mk_filt_hist(:,jj) = mkp1_plus;
            Pk = Pkp1_plus;
            Pk_filt_hist(:,:,jj)= Pkp1_plus;
            innovk_hist(:,jj) = innov_kp1;
            
            %%Compute and store NEES and NIS statistics:
            invPkp1 = inv(Pkp1_plus);
            invPyykp1 = inv(Pyykp1);
            NEESsshist(jj) = (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
            NISsshist(jj) = innov_kp1'*invPyykp1*innov_kp1;    
        end
    end
    NEESsamps(ss,:) = NEESsshist;
    NISsamps(ss,:) = NISsshist;
    
    %%sanity check:
    %Plot state estimation errors versus time
    f=figure(5),
    
    subplot(211)
    plot(t,xk_truehist(1,:)-mk_filt_hist(1,:),'Color',[0.5 0.5 0.5],'LineWidth',1), hold on
    plot(t,2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'m','LineWidth',1)
    plot(t,-2*sqrt(squeeze(Pk_filt_hist(1,1,:))'),'m','LineWidth',1)
    ylabel('$\xi$ (m)','Interpreter','latex')
    grid on, title('State Estimation Errors')
    subplot(212)
    plot(t,xk_truehist(2,:)-mk_filt_hist(2,:),'Color',[0.5 0.5 0.5],'LineWidth',1), hold on
    plot(t,2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'m','LineWidth',1)
    plot(t,-2*sqrt(squeeze(Pk_filt_hist(2,2,:))'),'m','LineWidth',1)
    grid on
    xlabel('Time (s)')
    ylabel('$\dot \xi$ (m)','Interpreter','latex')
end

