function [] = plotMeasurements(t, us, vs, lmks_visible, lmksToPlot)
%PLOTMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

% shapes = ['o','+','*','.','x','square',"diamond",'v','^','>',"Pentagram","Hexagram", "Horizontal line", "Vertical line"];
shapes = ['o','+','*','.','x','square',"diamond",'v','^','>'];

figure
hold on
for i= 1:length(lmksToPlot)
    lmk = lmksToPlot(i);
    u_plot = lmks_visible(lmk,:).*us(lmk,:);
    u_plot(u_plot==0) = nan;
    scatter(t,u_plot,shapes(mod(i-1,length(shapes))+1));
end
legend(split(int2str(lmksToPlot)))

figure
hold on
for i = 1:length(lmksToPlot)
    lmk = lmksToPlot(i);
    v_plot = lmks_visible(lmk,:).*vs(lmk,:);
    v_plot(v_plot==0) = nan;
    scatter(t,v_plot,shapes(mod(i-1,length(shapes))+1))
    
end

legend(split(int2str(lmksToPlot)))
end

