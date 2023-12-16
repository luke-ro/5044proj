function [Y,lmks_vis] = parseYtable(ytable)
%PARSEYTABLE Summary of this function goes here
%   Detailed explanation goes here

t_obs = unique(ytable(:,1));
lmks_vis = zeros(50,length(t_obs));
Y = NaN(100,length(t_obs));
for i = 1:length(t_obs)
    t = t_obs(i);
    
    mini_tab = ytable(ytable(:,1)==t,:);
    idxs = mini_tab(:,2);
    lmks_vis(idxs,i) = 1;
    
    Y(2*(idxs-1)+1,i) = mini_tab(:,3);
    Y(2*(idxs-1)+2,i) = mini_tab(:,4);
end
end

