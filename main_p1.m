% main_p1.m: part one main, simulate data
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
clc
close all

addpath('./functions/');
addpath('./Data/');
load('orbitdetermination-finalproj_data_2023_11_14.mat')

% Constants
% _N indicates a vector in the asteroid inertial frame
mu = 4.892*10^-9; % [km^3/s^2] Asteroid Gravitational Parameter
omegaA = 2*pi/(4.296057*60^2); % [s] Asteroid Rotation Rate
% omegaA = 4.296057*60^2; % [rad/s] Asteroid Orbital Period
rSA_N = [1.5*10^8; 0; 0]; % Inertial position of the sun wrt asteroid center
rSA = norm(rSA_N); % Distance from Sun to Asteroid
phi0 = 1*10^14; % [kg*km/s^2]
rho = 0.4; % coefficient of reflectivity
AreaMass = 1/62*10^-6; % Area-to-mass ratio

%% Problem 1:
% Simulate data with nonlinear dynamics and no noise

Dt = 600; % [s]
time_span = 0:Dt:259200; % [s]
params = [0,0]; %if we need extra constants for our func
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

r0_nom_N = [0,-1,0]';
r0_nom = norm(r0_nom_N);
v0_nom_N = [0,0,sqrt(mu/r0_nom)]';
X0_nom_N = [r0_nom_N; v0_nom_N];

ode_fun = @(t,X) dynamics(t,X,params);
[t,X_sim_N] = ode45(ode_fun,time_span,X0_nom_N,ode_options);
X_sim_N = X_sim_N';
t = t';

% plot inertial orbit
title = "True Trajectory in Inertial Frame";
plotOrbit(X_sim_N(1:3,:), title)

title = "States vs. Time";
plotStates(t,X_sim_N,title)

% DCMs
NC = R_CtoN;
NA = zeros(3,3,433);
X_sim_A = zeros(6,433);
AN = zeros(3,3,433);
for i = 1:length(t)
    theta = omegaA*t(i);
    NA(:,:,i) = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0 0 1];
    AN(:,:,i) = NA(:,:,i)';
    X_sim_A(:,i) = blkdiag(AN(:,:,i),AN(:,:,i))*X_sim_N(:,i);

%     CN(:,:,i) = NC(:,:,i)';
end
r_A = X_sim_A(1:3,:);

plotOrbit(r_A, "True Trajectory in Asteroid Frame")

%% Problem 2 jacobians
% derive jacobians for the dynamics


%% Problem 3 linearize
% linearize our dynamics about a given point


%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case
