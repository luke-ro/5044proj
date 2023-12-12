function [outputArg1,outputArg2] = calcNEES(n,inputArg2)
%CALCNEES Summary of this function goes here
%   Detailed explanation goes here
for i = 1:n
    % calculate initial state using m0 and P0:
    X0 = mvnrnd(x0,P0);


    %calculate true trajectory and measurements using non-linear dynamics
    X_true = 
    y_true = 

    % run the kalman filter 
    % X_est: (n x time) matrix of estimated states
    % P_yy:  (p x p X time) matrix of Pyy=0.5*H_kp1*P_minus_kp1*H_kp1'
    % P_plus: (n x n x time) matrix of covariance matrix after meas update
    % y_innov: (p x n x time) matrix of dynamics predicted y - meas y
    [X_est, P_yy, P_plus, y_innov] = kf();



    for j = 1:length(X_est)
        % NEES and NIS eqs from given code
        %     NEESsshist(jj) = (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
        %     NISsshist(jj) = innov_kp1'*invPyykp1*innov_kp1;  
        NEES_hist(j,i) = X_true(:,j) - X_est(:,j)'*(P_plus(:,:,j)^-1)*(X_true(:,j) - X_est(:,j);
        NIS_hist(j,i) = y_innov(:,j)'*P_yy(:,:,j)*y_innov(:,i);
    end
%     NEESsshist(jj) = (xk_truehist(:,jj) - mkp1_plus)'*invPkp1*(xk_truehist(:,jj) - mkp1_plus);
%     NISsshist(jj) = innov_kp1'*invPyykp1*innov_kp1;  
end

end

