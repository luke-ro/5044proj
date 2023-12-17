function [] = plotFilterResults(t_int, t_obs, X_sim, X_est, P_est, title, modifier)
%PLOTFILTERRESULTS Summary of this function goes here
%   Detailed explanation goes here
t_int = t_int/60^2;
t_obs = t_obs/60^2;

sgtitle(title, 'Interpreter', 'latex')
subplot(6,1,1)
hold on
plot(t_obs,X_est(1,:))
plot(t_obs,X_est(1,:)+2*sqrt(squeeze(P_est(1,1,:))'), '--')
plot(t_obs,X_est(1,:)-2*sqrt(squeeze(P_est(1,1,:))'), '--')
if length(X_sim) > 1
    plot(t_int, X_sim(1,:))
end
ylabel(strcat(modifier,'$r_{\hat{X}}$ [km]'), 'Interpreter', 'latex')
xlim([t_obs(2), t_obs(end)])

subplot(6,1,2)
hold on
plot(t_obs,X_est(2,:))
plot(t_obs,X_est(2,:)+2*sqrt(squeeze(P_est(2,2,:))'), '--')
plot(t_obs,X_est(2,:)-2*sqrt(squeeze(P_est(2,2,:))'), '--')
if length(X_sim) > 1
    plot(t_int, X_sim(2,:))
end
ylabel(strcat(modifier,'$r_{\hat{Y}}$ [km]'), 'Interpreter', 'latex')
xlim([t_obs(2), t_obs(end)])

subplot(6,1,3)
hold on
plot(t_obs,X_est(3,:))
plot(t_obs,X_est(3,:)+2*sqrt(squeeze(P_est(3,3,:))'), '--')
plot(t_obs,X_est(3,:)-2*sqrt(squeeze(P_est(3,3,:))'), '--')
if length(X_sim) > 1
    plot(t_int, X_sim(3,:))
end
ylabel(strcat(modifier,'$r_{\hat{Z}}$ [km]'), 'Interpreter', 'latex')
xlim([t_obs(2), t_obs(end)])

subplot(6,1,4)
hold on
plot(t_obs,X_est(4,:))
plot(t_obs,X_est(4,:)+2*sqrt(squeeze(P_est(4,4,:))'), '--')
plot(t_obs,X_est(4,:)-2*sqrt(squeeze(P_est(4,4,:))'), '--')
if length(X_sim) > 1
    plot(t_int, X_sim(4,:))
end
ylabel(strcat(modifier,'$v_{\hat{X}}$ [km/s]'), 'Interpreter', 'latex')
xlim([t_obs(2), t_obs(end)])

subplot(6,1,5)
hold on
if length(X_sim) > 1
    plot(t_int, X_sim(5,:))
end
plot(t_obs,X_est(5,:))
plot(t_obs,X_est(5,:)+2*sqrt(squeeze(P_est(5,5,:))'), '--')
plot(t_obs,X_est(5,:)-2*sqrt(squeeze(P_est(5,5,:))'), '--')
ylabel(strcat(modifier,'$v_{\hat{Y}}$ [km/s]'), 'Interpreter', 'latex')
xlim([t_obs(2), t_obs(end)])

subplot(6,1,6)
hold on
plot(t_obs,X_est(6,:))
plot(t_obs,X_est(6,:)+2*sqrt(squeeze(P_est(6,6,:))'), '--')
plot(t_obs,X_est(6,:)-2*sqrt(squeeze(P_est(6,6,:))'), '--')
if length(X_sim) > 1
    plot(t_int, X_sim(6,:))
end
ylabel(strcat(modifier,'$v_{\hat{Z}}$ [km/s]'), 'Interpreter', 'latex')
xlabel('time [hours]')
legend("$\delta x$", "$+2\sigma$","$+2\sigma$", "True $\delta x$","Interpreter", "Latex", "Location", "best")
xlim([t_obs(2), t_obs(end)])
end

