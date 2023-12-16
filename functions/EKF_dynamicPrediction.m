function [xhat_kplus1_minus, P_kplus1_minus] = EKF_dynamicPrediction(xhat_k_plus, P_k_plus, Q, gamma, const, steps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    ode_options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
    ode_fun = @(t,X) dynamics(t,X,const,zeros(length(xhat_k_plus),steps));
    [~,X_sim] = ode45(ode_fun,[0 const.Dt_int],xhat_k_plus,ode_options);
%     xhat_kplus1_minus = X_sim_N(:,1);
    xhat_kplus1_minus = X_sim(end,:)';

    %calcuate jacobian accoring to nominal traj
    Atilde = dyn_jacobian(xhat_kplus1_minus,const);
    
    %calcuate f using matrix expm
    Ahat = [Atilde zeros(6,1); zeros(1,7)];
    phim_hat = expm(Ahat*const.Dt_int);
    Ftilde = phim_hat(1:6, 1:6);
    Gtilde = phim_hat(1:6, 7);
    omega = const.Dt_int * gamma;

    P_kplus1_minus = Ftilde*P_k_plus*Ftilde' + omega*Q*omega';
    P_kplus1_minus = (P_kplus1_minus + P_kplus1_minus')/2;
end