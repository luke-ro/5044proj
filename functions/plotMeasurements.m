function [] = plotMeasurements(t, us, vs, lmks_visible, lmksToPlot, plot_title)
%PLOTMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here

t = t/60^2;

shapes = ['o','+','*','.','x','square',"diamond",'v','^','>','Pentagram',"Hexagram"];
% shapes = ['o','+','*','.','x','square',"diamond",'v','^','>'];

figure
subplot(2,1,1)
sgtitle(plot_title)
hold on
for i= 1:length(lmksToPlot)
    lmk = lmksToPlot(i);
    u_plot = lmks_visible(lmk,:).*us(lmk,:);
    u_plot(u_plot==0) = nan;
    scatter(t,u_plot,shapes(mod(i-1,length(shapes))+1));
end
legend(split(int2str(lmksToPlot)))
ylabel("u [pixels]")

subplot(2,1,2)
hold on
for i = 1:length(lmksToPlot)
    lmk = lmksToPlot(i);
    v_plot = lmks_visible(lmk,:).*vs(lmk,:);
    v_plot(v_plot==0) = nan;
    scatter(t,v_plot,shapes(mod(i-1,length(shapes))+1))
    
end

legend(split(int2str(lmksToPlot)))
xlabel("time [hours]")
ylabel("v [pixels]")
end

