function [X_sim_N, t, X_simObs_N, t_obs] = simNLdynamics(w_tilde, X0, const)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% time_span = 0:const.Dt_int:const.tf_int; % [s]
% time_span = 0:const.Dt_obs:const.tf_obs; % [s]
time_span = 0:const.Dt_int:const.tf_obs; % [s]
% X0_nom_N = [const.r0_nom_N; const.v0_nom_N];
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
ode_fun = @(t,X) dynamics(t,X,const,w_tilde);
[t,X_sim_N] = ode45(ode_fun,time_span,X0,ode_options);
X_sim_N = X_sim_N';
t = t';

X_simObs_N = X_sim_N(:,1:10:length(t));
t_obs = t(1:10:length(t));



end

