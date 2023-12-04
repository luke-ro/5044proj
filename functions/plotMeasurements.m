function [] = plotMeasurements(t, us, vs, lmks_visible, lmksToPlot)
%PLOTMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

shapes = ['o','+','*','.','x','square',"diamond",'v','^','>'];

figure
hold on
for lmk = lmksToPlot
    u_plot = lmks_visible(lmk,:).*us(lmk,:);
    u_plot(u_plot==0) = nan;
    scatter(t,u_plot,shapes(lmk))
end
figure
hold on
for lmk = lmksToPlot
    v_plot = lmks_visible(lmk,:).*vs(lmk,:);
    v_plot(v_plot==0) = nan;
    scatter(t,v_plot,shapes(lmk))
end
end

