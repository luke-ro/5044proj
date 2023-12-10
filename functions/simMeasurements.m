function [us, vs, lmks_visible] = simMeasurements(t, X_sim_N, NC, pos_lmks_A, const)
%SIMMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

num_LMKs = 50;
% n = 6;
% npts = length(t);

% NA = zeros(3,3,433);
X_sim_A = zeros(6,433);
% AN = zeros(3,3,433);
% pos_lmks_C = zeros(3,50,433);
% pos_lmks_N = zeros(3,50,433);

lmks_in_FOV = zeros(num_LMKs,length(t));
lmks_in_front = zeros(num_LMKs,length(t));

us = zeros(num_LMKs,length(t));
vs = zeros(num_LMKs,length(t));

for i = 1:length(t)
    theta = const.omegaA*t(i);
    NA = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0 0 1];
    AN = NA';
%     X_sim_A(:,i) = blkdiag(AN(:,:,i),AN(:,:,i))*X_sim_N(:,i);
    r_A = AN*X_sim_N(1:3,i);
    v_A = AN*X_sim_N(4:6,i);
    X_sim_A(:,i) = [r_A;v_A];
    
    CN = NC(:,:,i)';
    AC = AN*NC(:,:,i);
    
    % get a logical vec of true where landmarks are in "front" of benu
    lmks_in_front(:,i) = getLMsInFront(CN,NA,pos_lmks_A);

    % get the pixel coords of all landmarks
    [us(:,i),vs(:,i)] = uv_func(r_A, pos_lmks_A, AC, const.u0, const.v0);
    
    % add noise
%     R = diag([const.sig_uv^2, const.sig_uv^2]);
    v_u = mvnrnd(0, const.sig_uv^2);
    v_v = mvnrnd(0, const.sig_uv^2);
    us(:,i) = us(:,i) + v_u;
    vs(:,i) = vs(:,i) + v_v;
    
    % get a logical vec of true where LMs are within FOV
    lmks_in_FOV(:,i) = getLMsInFOV(pos_lmks_A, r_A, AC(:,3), us(:,i)', vs(:,i)',const.uv_max,const.uv_max);

%     pos_lmks_C(:,:,i) = AC'*pos_lmks_A;
%     pos_lmks_N(:,:,i) = NA*pos_lmks_A;
end

% get logical vec vector true where lmks are visible to satellite
lmks_visible = lmks_in_FOV & lmks_in_front;

end

