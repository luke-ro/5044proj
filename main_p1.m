% main_p1.m: part one main, simulate data
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
clc

addpath('./functions/');
addpath('./Data/');
load('orbitdetermination-finalproj_data_2023_11_14.mat')

% Constants
% _N indicates a vector in the asteroid inertial frame
mu = 4.892*10^-9; %[km^3/s^2] Asteroid Gravitational Parameter
rSA_N = [1.5*10^8; 0; 0]; % Inertial position of the sun wrt asteroid center
rSA = norm(rSA_N); % Distance from Sun to Asteroid
phi0 = 1*10^14; % [kg*km/s^2]
rho = 0.4; % coefficient of reflectivity
AreaMass = 1/62*10^-6; % Area-to-mass ratio

%% Problem 1:
% Simulate data with nonlinear dynamics and no noise

Dt = 600; % [s]
time_span = [0,100000]; % [s]
params = [0,0]; %if we need extra constants for our func
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

r0_nom_N = [0,-1,0]';
r0_nom = norm(r0_nom_N);
v0_nom_N = [0,0,sqrt(mu/r0_nom)]';
x0_nom_N = [r0_nom_N; v0_nom_N];

ode_fun = @(t,x) dynamics(t,x,params);
[t,X_sim] = ode45(ode_fun,time_span,x0_nom_N,ode_options);

% plot inertial position, earth, and max altitude location
figure
plot3(X_sim(:,1), X_sim(:,2), X_sim(:,3))
hold on
% surf(surfX, surfY, surfZ)
axis equal
title("Satellite Orbit around Asteroid (Inertial Frame)")
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off
%% Problem 2 jacobians
% derive jacobians for the dynamics


%% Problem 3 linearize
% linearize our dynamics about a given point


%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case
