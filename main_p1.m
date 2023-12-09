% main_p1.m: part one main, simulate data
% ASEN 5044 Fall 23 Project
% Mikaela Felix, Amanda Marlow, Luke Roberson

clear
close all

addpath('./functions/');
addpath('./Data/');
load('orbitdetermination-finalproj_data_2023_11_14.mat')

% Constants
% _N indicates a vector in the asteroid inertial frame
const.mu = 4.892*10^-9; % [km^3/s^2] Asteroid Gravitational Parameter
const.omegaA = 2*pi/(4.296057*60^2); % [s] Asteroid Rotation Rate
const.rSA_N = [1.5*10^8; 0; 0]; % Inertial position of the sun wrt asteroid center
const.rSA = norm(const.rSA_N); % Distance from Sun to Asteroid
const.phi0 = 1*10^14; % [kg*km/s^2]
const.rho = 0.4; % coefficient of reflectivity
const.AreaMass = 1/62*10^-6; % Area-to-mass ratio
const.f = 2089.7959; % [pixels]
const.u0 = 512; % [pixels]
const.v0 = 512; % [pixels]
const.uv_min = 0; % [pixels]
const.uv_max = 1024; % [pixels]

num_LMKs = 50;

%% Problem 1:
% Simulate data with nonlinear dynamics and no noise

Dt = 600; % [s]
time_span = 0:Dt:259200; % [s]
params = [0,0]; %if we need extra constants for our func
ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

r0_nom_N = [0,-1,0]';
r0_nom = norm(r0_nom_N);
v0_nom_N = [0,0,sqrt(const.mu/r0_nom)]';
X0_nom_N = [r0_nom_N; v0_nom_N];

ode_fun = @(t,X) dynamics(t,X,params);
[t,X_nom_N] = ode45(ode_fun,time_span,X0_nom_N,ode_options);
X_nom_N = X_nom_N';
t = t';

% plot inertial orbit
title = "True Trajectory in Inertial Frame";
plotOrbit(X_nom_N(1:3,:), title)

figure
title = "States vs. Time, Full Nonlinear Dynamics Simulation";
plotStates(t,X_nom_N,title,"")

% DCMs
NC = R_CtoN;
NA = zeros(3,3,433);
X_sim_A = zeros(6,433);
AN = zeros(3,3,433);
pos_lmks_C = zeros(3,50,433);
pos_lmks_N = zeros(3,50,433);
uv = zeros(2,50,433);


lmks_in_FOV = zeros(num_LMKs,length(t));
lmks_in_front = zeros(num_LMKs,length(t));

us = zeros(num_LMKs,length(t));
vs = zeros(num_LMKs,length(t));
for i = 1:length(t)
    theta = const.omegaA*t(i);
    NA(:,:,i) = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0 0 1];
    AN(:,:,i) = NA(:,:,i)';
%     X_sim_A(:,i) = blkdiag(AN(:,:,i),AN(:,:,i))*X_sim_N(:,i);
    r_A = AN(:,:,i)*X_nom_N(1:3,i);
    v_A = AN(:,:,i)*X_nom_N(4:6,i);
    X_sim_A(:,i) = [r_A;v_A];
    
    CN(:,:,i) = NC(:,:,i)';
    AC(:,:,i) = AN(:,:,i)*NC(:,:,i);
    
    % get a logical vec of true where landmarks are in "front" of benu
    lmks_in_front(:,i) = getLMsInFront(CN(:,:,i),NA(:,:,i),pos_lmks_A);

    % get the pixel coords of all landmarks
    [us(:,i),vs(:,i)] = uv_func(r_A, pos_lmks_A, AC(:,:,i),const.u0,const.v0);
    
    % get a logical vec of true where LMs are within FOV
    lmks_in_FOV(:,i) = getLMsInFOV(pos_lmks_A, r_A, AC(:,3,i),us(:,i)',vs(:,i)',const.uv_max,const.uv_max);

    pos_lmks_C(:,:,i) = AC(:,:,i)'*pos_lmks_A;
    pos_lmks_N(:,:,i) = NA(:,:,i)*pos_lmks_A;
end

% get logical vec vector true where lmks are visible to satellite
lmks_visible = lmks_in_FOV & lmks_in_front;

% get position from simulated points
r_A = X_sim_A(1:3,:);

plotOrbit(r_A, "True Trajectory in Asteroid Frame")
scatter3(pos_lmks_A(1,:), pos_lmks_A(2,:), pos_lmks_A(3,:), '.')

plotMeasurements(t, us, vs, lmks_visible, 1:10, "Full Nonlinear Measurement Simulation")

%% Problem 2 jacobians
% derive jacobians for the dynamics

%% Problem 3 linearize
% linearize our dynamics about a given point

% B and D matrices
B = zeros(6,1);
% B(4:6) = -const.phi0/const.rSA^3*(1+4/9*const.rho)*const.AreaMass * const.rSA_N; % <- solution doesn't include this
D = zeros(2,2);

% initialize A and C matrices
bigA = zeros(6,6,length(time_span));
bigF = zeros(6,6,length(time_span));
bigG = zeros(6,1,length(time_span));
bigC = cell(1,length(time_span));
bigCsym = cell(1,length(time_span));

% C = zeros(2, 6, time_span);
z_1by7 = zeros(1,7);
X0_delta = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';

% X0_delta = [0 0 0 1e7  0 0]';
X_delta_N = [X0_delta zeros(6,length(time_span)-1)];
X_DT_total = [X0_nom_N+X0_delta zeros(6,length(time_span)-1)];
X_DT_nom = [X0_nom_N zeros(6,length(time_span) - 1)];
Y_delta_N = cell(1,length(time_span));
u_delta_N = NaN(length(pos_lmks_A),length(time_span));
v_delta_N = NaN(length(pos_lmks_A),length(time_span));

sym_jacobian_H = makeSymH(const);

for i = 1:length(time_span)-1
    %calcuate jacobian accoring to nominal traj
    A = dyn_jacobian(X_nom_N(:,i),const);
    
    %calcuate f using matrix expm
    Ahat = [A B; z_1by7];
    phim_hat = expm(Ahat*Dt);
    F = phim_hat(1:6, 1:6);
    G = phim_hat(1:6, 7);
    
%     F = eye(6) +A*Dt;

    %propogate delta state
    X_delta_N(:,i+1) = F * X_delta_N(:,i) + G;
    
    % cacluate C
    num_landmarks = sum(lmks_visible(:,i));
    lmk_idxs = find(lmks_visible(:,i));
    num_measurements = num_landmarks*2;% (number of landmarks in view)*2 measurements
    C = zeros(num_measurements,6); 
    C_sym = zeros(num_measurements,6); 
    for j = 1:2:num_measurements
        C(j:j+1,:) = dyn_jacobian_H(X_nom_N(:,i), pos_lmks_N(:,lmk_idxs((j+1)/2),i), NC(:,:,i));
        C_sym(j:j+1,:) = sym_jacobian_H(X_nom_N(:,i), pos_lmks_N(:,lmk_idxs((j+1)/2),i), NC(:,:,i));      
    end
    
    Y_delta_N(i) = {C_sym*X_delta_N(:,i)};
    u_delta_N(lmk_idxs,i) = Y_delta_N{i}(1:2:end);
    v_delta_N(lmk_idxs,i) = Y_delta_N{i}(2:2:end);
    
    bigA(:,:,i) = A;
    bigF(:,:,i) = F;
    bigG(:,:,i) = G;
    bigC(i) = {C};
    bigCsym(i) = {C_sym};
end
figure;hold on;
plotStates(t,X_delta_N,"State Deviations vs. Time, Linearized Dynamics Simulation","$$\delta$$")


figure;hold on;
scatter(t/3600,u_delta_N(1,:))
%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case


save("data/P1_vars.mat","X_nom_N","X_delta_N","Y_delta_N","bigA","bigF","bigC","pos_lmks_N","lmks_visible");


