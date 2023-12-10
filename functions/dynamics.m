function [dX] = dynamics(t, X, const, w_tilde)
% This is the function ode45 uses to integrate the nonlinear dynamics

% Constants
% _N indicates a vector in the asteroid inertial frame
% mu = 4.892*10^-9; %[km^3/s^2] Asteroid Gravitational Parameter
% rSA_N = [1.5*10^8; 0; 0]; % Inertial position of the sun wrt asteroid center
% rSA = norm(rSA_N); % Distance from Sun to Asteroid
% % rSA_uv_N = rsa_N/
% phi0 = 1*10^14; % [kg*km/s^2]
% rho = 0.4; % coefficient of reflectivity
% AreaMass = 1/62*10^-6; % Area-to-mass ratio

% Extract states
r_N = X(1:3); % inertial position
r = norm(r_N); % radial distance
v_N = X(4:6); % inertial velocity

a2B = -const.mu/(r^3)*r_N; % two body acceleration
aSRP = -const.phi0/const.rSA^3*(1+4/9*const.rho)*const.AreaMass * const.rSA_N; % solar radiation pressure
a_N = a2B + aSRP; % total acceleration

dX = [v_N; a_N];

% add noise
if t > 0
%     dX = dX + w_tilde(:,ceil(t/const.Dt_int));
    dX = dX + w_tilde(:,ceil(t/const.Dt_obs));
end

end

