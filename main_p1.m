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

    % the next three lines should replace the commented for loop below
    
    % get a logical vec of true where landmarks are in "front" of benu
    lmk_in_front_new = getLMsInFront(CN(:,:,i),NA(:,:,i),pos_lmks_A);

    %get the pixel coords of all landmarks
    [us,vs] = uv_func(r_A, pos_lmks_A, NC(:,:,i),u0,v0);
    
    % get a logical vec of true where LMs are within FOV
    lmks_in_FOV_new = getLMsInFOV(pos_lmks_A,r_A,NC(:,3,i),us,vs,uv_max,uv_max);

%     for j = 1:length(pos_lmks_A)
%         pos_lmks_C(:,j,i) = CN(:,:,i)*NA(:,:,i)*pos_lmks_A(:,j);
%         lmk_in_front(i,j) = pos_lmks_C(3,j,i) < 0;
% 
%         [u,v] = uv_func(r_A, pos_lmks_A(:,j), NC(:,:,i),u0,v0);
%         lmk_in_FOV = (0<=u & u<=uv_max) & (0<=v & v<=uv_max) & (dot((pos_lmks_A(:,j) - r_A),NC(:,3,i)) > 0);
%         if pos_lmks_C(3,j,i) < 0
%             uv(:,j,i) = [u;v];
%         end
%     end
end
r_A = X_sim_A(1:3,:);

plotOrbit(r_A, "True Trajectory in Asteroid Frame")
scatter3(pos_lmks_A(1,:), pos_lmks_A(2,:), pos_lmks_A(3,:), '.')

%% Problem 2 jacobians
% derive jacobians for the dynamics


%% Problem 3 linearize
% linearize our dynamics about a given point

% B and D matrices
B = [0; 0; 0; -phi0/rSA^3*(1+4/9*rho)*AreaMass * rSA_N];
D = zeros(2,2);

% initialize A and C matrices
bigA = zeros(6,6,length(time_span));
bigF = zeros(6,6,length(time_span));
bigG = zeros(6,1,length(time_span));

% C = zeros(2, 6, time_span);
z_1by7 = zeros(1,7);
X0_delta = [1e-5 1e-5 1e-5 1e-7 1e-7 1e-7]';
X_delta = [X0_delta zeros(6,length(time_span)-1)];
X_DT_total = [X0_nom_N+X0_delta zeros(6,length(time_span)-1)];
X_DT_nom = [X0_nom_N zeros(6,length(time_span) - 1)];

for i = 1:length(time_span)
    A = dyn_jacobian(X_sim_N(:,i));
    Ahat = [A B; z_1by7];
    phim_hat = expm(Ahat*Dt);
    F = phim_hat(1:6, 1:6);
    G = phim_hat(1:6, 7);
    
    X_DT_nom(:,i+1) = F * X_DT_nom(:,i) + G;
    X_DT_total(:,i+1) = F * X_DT_total(:,i);

%     X_DT_total(:,i+1) = F * (X_sim_N(:,i)+X_delta(:,i));
    X_delta(:,i+1) = X_DT_total(:,i+1) - X_DT_nom(:,i);

<<<<<<< Updated upstream
    bigA(:,:,i) = A;
    bigF(:,:,i) = F;
    bigG(:,:,i) = G;

    for j = 1:length(pos_lmks_A)
        C = dyn_jacobian_H(X_sim_N(:,i), pos_lmks_C(:,j,i), R_CtoN(:,:,i));
    end
=======
    %propogate delta state
    X_delta_N(:,i+1) = F * X_delta_N(:,i) + G;
    
    % cacluate C
    num_landmarks = sum(lmks_visible(:,i));
    lmk_idxs = find(lmks_visible(:,i));
    num_measurements = num_landmarks*2;% (number of landmarks in view)*2 measurements
    C = zeros(num_measurements,6); 
    C_sym = zeros(num_measurements,6); 
    for j = 1:2:num_measurements
        C(j:j+1,:) = dyn_jacobian_H(X_nom_N(:,i), pos_lmks_N(:,lmk_idxs((j+1)/2),i), NC(:,:,i), const);
%         C_sym(j:j+1,:) = sym_jacobian_H(X_nom_N(:,i), pos_lmks_N(:,lmk_idxs((j+1)/2),i), NC(:,:,i));      
    end
    
    Y_delta_N(i) = {C*X_delta_N(:,i)};
    u_delta_N(lmk_idxs,i) = Y_delta_N{i}(1:2:end);
    v_delta_N(lmk_idxs,i) = Y_delta_N{i}(2:2:end);
    
    bigA(:,:,i) = A;
    bigF(:,:,i) = F;
    bigG(:,:,i) = G;
    bigC(i) = {C};
%     bigCsym(i) = {C_sym};
>>>>>>> Stashed changes
end
X_DT_nom = X_DT_nom(:,1:end-1);
X_DT_total = X_DT_total(:,1:end-1);

% X_DT_nom = [X0_nom_N+X0_delta zeros(6,length(time_span) - 1)];
% for i = 1:length(time_span)
%     A = dyn_jacobian(X_DT_nom(:,i));
%     Ahat = [A B; z_1by7];
%     phim_hat = expm(Ahat*Dt);
%     F = phim_hat(1:6, 1:6);
%     G = phim_hat(1:6, 7);
%     X_delta(:,i+1) = F * X_delta(:,i);
%     X_DT_nom(:,i+1) = F * X_DT_nom(:,i);
%     
% 
%     bigA(:,:,i) = A;
%     bigF(:,:,i) = F;
%     bigG(:,:,i) = G;
% 
%     for j = 1:length(pos_lmks_A)
%         C = dyn_jacobian_H(X_sim_N(:,i), pos_lmks_C(:,j,i), R_CtoN(:,:,i));
%     end
% end
% X_DT_nom = X_DT_nom(:,2:end);
X_delta = X_delta(:,2:end);

% figure
% hold on
% sgtitle("Nominal State over Time")
% 
% subplot(6,1,1); plot(t,X_DT_nom(1,:))
% 
% subplot(6,1,2); plot(t,X_DT_nom(2,:))
% 
% subplot(6,1,3); plot(t,X_DT_nom(3,:))
% 
% subplot(6,1,4); plot(t,X_DT_nom(4,:))
% 
% subplot(6,1,5); plot(t,X_DT_nom(5,:))
% 
% subplot(6,1,6); plot(t,X_DT_nom(6,:))
% 
% hold off

% plotStates(t,X_DT_nom,"Nominal State over Time")




%% Problem 4 Compare nonlinear to linear
% simulate linearized dynamics and compare to the nonlinear case

figure
hold on
sgtitle("Nominal State over Time")

subplot(6,1,1)
hold on
plot(t,X_DT_nom(1,:))
plot(t,X_sim_N(1,:))
ylabel("x")
hold off

subplot(6,1,2); 
hold on
plot(t,X_DT_nom(2,:))
plot(t,X_sim_N(2,:))
ylabel("y")
hold off

subplot(6,1,3); 
hold on
plot(t,X_DT_nom(3,:))
plot(t,X_sim_N(3,:))
ylabel("z")
hold off

subplot(6,1,4); 
hold on
plot(t,X_DT_nom(4,:))
plot(t,X_sim_N(4,:))
ylabel("xdot")
hold off

subplot(6,1,5);
hold on
plot(t,X_DT_nom(5,:))
plot(t,X_sim_N(5,:))
ylabel("ydot")
hold off

subplot(6,1,6); 
hold on
plot(t,X_DT_nom(6,:))
plot(t,X_sim_N(6,:))
ylabel("zdot")
hold off

hold off

