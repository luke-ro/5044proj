% main_p1.m: part one main, simulate data
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson
addpath('./functions/');

%% Problem 1:
% Simulate data with nonlinear dynamics and no noise

Dt = 600; % [s]
time_span = [0,100000]; % [s]
params = [0,0]; %if we need extra constants for our func
ode_options = odeset('RelTol',1e-12);

x0 = [0,0,0,0,0,0];

ode_fun = @(t,x) dynamics(t,x,params);
[t,X_sim] = ode45(ode_fun,time_span,x0,ode_options);

%% Problem 2 jacobians
% derive jacobians for the dynamics


%% Problem 3 linearize
% linearize our dynamics about a given point


%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case
