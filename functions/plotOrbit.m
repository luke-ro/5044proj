function [] = plotOrbit(r_V, figTitle)
%PLOTORBIT Summary of this function goes here
%   Detailed explanation goes here

figure
plot3(r_V(1,:), r_V(2,:), r_V(3,:), 'LineWidth', 1.5)
hold on
scatter3(0,0,0, 'filled')
axis equal
title(figTitle)
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
grid on

end

