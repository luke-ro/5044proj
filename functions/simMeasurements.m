function [us, vs, lmks_visible] = simMeasurements(t, X_sim_N, NC, pos_lmks_A, const)
%SIMMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

num_LMKs = 50;
% n = 6;
npts = length(t)-1;

X_sim_A = zeros(6,npts);

lmks_in_FOV = zeros(num_LMKs,npts);
lmks_in_front = zeros(num_LMKs,npts);

us = zeros(num_LMKs,npts);
vs = zeros(num_LMKs,npts);

for i = 2:length(t)
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
    lmks_in_front(:,i-1) = getLMsInFront(CN,NA,pos_lmks_A);

    % get the pixel coords of all landmarks
    [us(:,i-1),vs(:,i-1)] = uv_func(r_A, pos_lmks_A, AC, const.u0, const.v0);
    lmks_in_FOV1 = getLMsInFOV(pos_lmks_A, r_A, AC(:,3), us(:,i-1)', vs(:,i-1)',const.uv_max,const.uv_max);
    
    % add noise
%     R = diag([const.sig_uv^2, const.sig_uv^2]);
    v_u = mvnrnd(zeros(num_LMKs,1), const.sig_uv^2);
    v_v = mvnrnd(zeros(num_LMKs,1), const.sig_uv^2);
    
    us(:,i-1) = us(:,i-1) + v_u;
    vs(:,i-1) = vs(:,i-1) + v_v;
    
    % get a logical vec of true where LMs are within FOV
    lmks_in_FOV2 = getLMsInFOV(pos_lmks_A, r_A, AC(:,3), us(:,i-1)', vs(:,i-1)',const.uv_max,const.uv_max);

%     pos_lmks_C(:,:,i) = AC'*pos_lmks_A;
%     pos_lmks_N(:,:,i) = NA*pos_lmks_A;

    lmks_in_FOV(:,i-1) = lmks_in_FOV1 & lmks_in_FOV2;
end

% get logical vec vector true where lmks are visible to satellite
lmks_visible = lmks_in_FOV & lmks_in_front;

end

