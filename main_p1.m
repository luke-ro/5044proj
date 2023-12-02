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
f = 2089.8959; % [pixels]
u0 = 512; % [pixels]
v0 = 512; % [pixels]
uv_min = 0; % [pixels]
uv_max = 1024; % [pixels]

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
pos_lmks_C = zeros(3,50,433);
uv = zeros(2,50,433);
for i = 1:length(t)
    theta = omegaA*t(i);
    NA(:,:,i) = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0 0 1];
    AN(:,:,i) = NA(:,:,i)';
%     X_sim_A(:,i) = blkdiag(AN(:,:,i),AN(:,:,i))*X_sim_N(:,i);
    r_A = AN(:,:,i)*X_sim_N(1:3,i);
    v_A = AN(:,:,i)*X_sim_N(4:6,i);
    X_sim_A(:,i) = [r_A;v_A];
    
    CN(:,:,i) = NC(:,:,i)';

    for j = 1:length(pos_lmks_A)
        pos_lmks_C(:,j,i) = CN(:,:,i)*NA(:,:,i)*pos_lmks_A(:,j);
        lmk_in_front(i,j) = pos_lmks_C(3,j,i) < 0;

        [u,v] = uv_func(r_A, pos_lmks_A(:,j), NC(:,:,i),u0,v0);
        lmk_in_FOV = (0<=u & u<=uv_max) & (0<=v & v<=uv_max) & (dot((pos_lmks_A(:,j) - r_A),NC(:,3,i)) > 0);
        if pos_lmks_C(3,j,i) < 0
            uv(:,j,i) = [u;v];
        end
    end
end
r_A = X_sim_A(1:3,:);

plotOrbit(r_A, "True Trajectory in Asteroid Frame")
scatter3(pos_lmks_A(1,:), pos_lmks_A(2,:), pos_lmks_A(3,:), '.')

%% Problem 2 jacobians
% derive jacobians for the dynamics


%% Problem 3 linearize
% linearize our dynamics about a given point


%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case
