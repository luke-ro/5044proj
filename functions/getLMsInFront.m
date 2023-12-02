function [lmk_in_front] = getLMsInFront(CN,NA,lms)
%GETLMSINFRONT Summary of this function goes here
%   Detailed explanation goes here
CA = CN*NA;
pos_lmks_C = CA*lms;
lmk_in_front = pos_lmks_C(3,:) < 0;
end

