function [stacked] = stackUsVs(us,vs)
%STACKUSVS Summary of this function goes here
%   Detailed explanation goes here
stacked = zeros(size(us,1)*2,size(us,2));
for k = 1:size(us,2)
    for i = 1:size(us,1)
        rows = (2*(i-1)+1):(2*(i-1)+2);
        stacked(rows,k) = [us(i);vs(i)];
    end
end

