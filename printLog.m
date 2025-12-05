% printLog.m
% -------------------------------------------------------------------------
% Utility script to load a saved simulation log (log.mat) and visualize:
%   1) Hitch reference (desired) path vs. actual hitch path in XY-plane
%   2) Control inputs over time: drive force F_drive and yaw torque tau_r
%
% Expected structure inside log.mat:
%   log.t    : [N x 1] time vector
%   log.y    : [N x 2] actual hitch position [x_h, y_h]
%   log.yd   : [N x 2] desired hitch position [x_h, y_h]
%   log.u    : [N x 2] control inputs [F_drive, tau_r]
%
% Make sure that 'log.mat' is in the current MATLAB folder or on the path
% before running this script.
% -------------------------------------------------------------------------

% Load log from file
loadinglog = load('log.mat');
log        = loadinglog.log;

% =========================================================
% 1) Hitch actual vs desired path
% =========================================================
figure('Color','w'); hold on; grid on; axis equal;
plot(log.yd(:,1), log.yd(:,2), 'r--', 'LineWidth', 1.5);  % desired hitch path
plot(log.y(:,1),  log.y(:,2),  'b-',  'LineWidth', 1.5);  % actual hitch path
xlabel('x_h'); ylabel('y_h');
legend('desired hitch','actual hitch','Location','Best');
title('Hitch path: actual vs desired');

% =========================================================
% 2) Control inputs over time
% =========================================================
figure;
subplot(2,1,1);
plot(log.t, log.u(:,1));
ylabel('F\_drive'); grid on;
title('Control inputs');

subplot(2,1,2);
plot(log.t, log.u(:,2));
ylabel('\tau_r'); xlabel('t'); grid on;
