function [X_sim_N, t] = simNLdynamics(w_tilde, const)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

time_span = 0:const.Dt_int:const.tf_int; % [s]
X0_nom_N = [const.r0_nom_N; const.v0_nom_N];
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
ode_fun = @(t,X) dynamics(t,X,const,w_tilde);
[t,X_sim_N] = ode45(ode_fun,time_span,X0_nom_N,ode_options);
X_sim_N = X_sim_N';
t = t';
end

