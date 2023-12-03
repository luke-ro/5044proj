function [lmk_in_FOV] = getLMsInFOV(lms,r,kc,us,vs,umax,vmax)
%GETLMSINFOV Summary of this function goes here
%   Detailed explanation goes here
kc_vec = repmat(kc,1,size(lms,2));
lmk_in_FOV = (0<=us & us<=umax) & (0<=vs & vs<=vmax) & (dot((lms - r),kc_vec,1) > 0);
end

