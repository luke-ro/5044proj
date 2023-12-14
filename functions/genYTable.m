function [tab] = genYTable(t, us, vs, vis)
%GENYTABLE Summary of this function goes here
%   Detailed explanation goes here
n = length(t);
tab = zeros(n*size(us,1),4);
row = 1;
for i = 1:n
    n_lmks_vis = sum(vis(:,i));
    cur_lmk_idxs = find(vis(:,i));
    row_end = row + n_lmks_vis-1;
    mini_table = [t(i)*ones(n_lmks_vis,1), cur_lmk_idxs, us(cur_lmk_idxs,i), vs(cur_lmk_idxs,i)];
    tab(row:row_end,:) = mini_table;
    row = row_end+1;
end

tab = tab(1:row-1,:);

