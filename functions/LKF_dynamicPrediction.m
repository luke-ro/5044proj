function [delta_x, P] = LKF_dynamicPrediction(delta_x0, P0, F, Q, OMEGA)
%LKFSTEP Single dynamic prediction step of Linearized Kalman Filter

delta_x =  F*delta_x0;
% P = F*P0*F' + OMEGA*Q*OMEGA'; % double check that we don't need gamma
P = F*P0*F' + Q; % double check that we don't need gamma


% ensure it stays positive definite
P = 1/2*(P + P');

end

