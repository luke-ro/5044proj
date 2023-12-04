function [] = plotStates(t,X,title)
%PLOTSTATES Summary of this function goes here
%   Detailed explanation goes here

figure
subplot(6,1,1)
sgtitle(title, 'Interpreter', 'latex')
plot(t, X(1,:));
hold on
ylabel('$r_{\hat{X}}$ [m]', 'Interpreter', 'latex')
hold off

subplot(6,1,2)
plot(t, X(2,:));
hold on
ylabel('$r_{\hat{Y}}$ [m]', 'Interpreter', 'latex')
hold off

subplot(6,1,3)
plot(t, X(3,:));
hold on
ylabel('$r_{\hat{Z}}$ [m]', 'Interpreter', 'latex')
xlabel('time [s]')
hold off

subplot(6,1,4)
plot(t, X(4,:));
hold on
ylabel('$v_{\hat{X}}$ [m/s]', 'Interpreter', 'latex')
hold off

subplot(6,1,5)
plot(t, X(5,:));
hold on
ylabel('$v_{\hat{Y}}$ [m/s]', 'Interpreter', 'latex')
hold off

subplot(6,1,6)
plot(t, X(6,:));
hold on
ylabel('$v_{\hat{Z}}$ [m/s]', 'Interpreter', 'latex')
xlabel('time [s]')
hold off

end

