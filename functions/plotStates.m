function [] = plotStates(t,X,title,modifier)
%PLOTSTATES Summary of this function goes here
%   Detailed explanation goes here

t = t/60^2;


subplot(6,1,1)
sgtitle(title, 'Interpreter', 'latex')
plot(t, X(1,:));
hold on
ylabel(strcat(modifier,'$r_{\hat{X}}$ [km]'), 'Interpreter', 'latex')
% hold off

subplot(6,1,2)
plot(t, X(2,:));
hold on
ylabel(strcat(modifier,'$r_{\hat{Y}}$ [km]'), 'Interpreter', 'latex')
% hold off

subplot(6,1,3)
plot(t, X(3,:));
hold on
ylabel(strcat(modifier,'$r_{\hat{Z}}$ [km]'), 'Interpreter', 'latex')
xlabel('time [hours]')
% hold off

subplot(6,1,4)
plot(t, X(4,:));
hold on
ylabel(strcat(modifier,'$v_{\hat{X}}$ [km/s]'), 'Interpreter', 'latex')
% hold off

subplot(6,1,5)
plot(t, X(5,:));
hold on
ylabel(strcat(modifier,'$v_{\hat{Y}}$ [km/s]'), 'Interpreter', 'latex')
% hold off

subplot(6,1,6)
plot(t, X(6,:));
hold on
ylabel(strcat(modifier,'$v_{\hat{Z}}$ [km/s]'), 'Interpreter', 'latex')
xlabel('time [hours]')
% hold off

end

