%   log    : struct with fields
%              .t  [N x 1]
%              .y  [N x 2]  actual trailer com position [x_t, y_t]
%              .X  [N x 12] full state (can be used for hitch, etc.)
%              .u  [N x 2]  control inputs [F_drive, tau_r]
%   ref_to : struct with fields
%              .p_hitch [2 x N] TO hitch reference path
%              .p_tr    [2 x N] TO trailer reference path
%   note: here y follows the trailer com, not the hitch



% Load from file
data = load('log.mat');

if ~isfield(data,'log')
    error('log.mat does not contain variable ''log''.');
end
if ~isfield(data,'ref_to')
    error('log.mat does not contain variable ''ref_to''.');
end

log    = data.log;
ref_to = data.ref_to;

% transpose reference paths: 2xN -> N x 2
p_hitch_ref = ref_to.p_r.';   % [N x 2]
p_tr_ref    = ref_to.p_tr.';      % [N x 2]

% =========================================================
% 1) trailer: TO reference path vs actual trailer com path
% =========================================================
figure('Color','w'); hold on; grid on; axis equal;
plot(p_tr_ref(:,1), p_tr_ref(:,2), 'r--', 'LineWidth', 1.5);   % TO ref (trailer)
plot(log.y(:,1),    log.y(:,2),    'b-',  'LineWidth', 1.5);   % actual trailer com
xlabel('x_t');
ylabel('y_t');
legend('TO reference path (trailer)','actual trailer com','Location','Best');
title('trailer com path: TO reference vs actual');

%{
% =========================================================
% 1b) Trailer: TO reference path vs actual trailer path  (from log.X)
% =========================================================
if isfield(log,'X') && size(log.X,2) >= 8
    xt = log.X(:,7);   % trailer x (根据你的状态定义 [xr,yr,thetar,vxr,vyr,wr,xt,yt,...])
    yt = log.X(:,8);   % trailer y
    figure('Color','w'); hold on; grid on; axis equal;
    plot(p_tr_ref(:,1), p_tr_ref(:,2), 'r--', 'LineWidth', 1.5);   % TO ref (trailer)
    plot(xt, yt, 'b-', 'LineWidth', 1.5);                          % actual trailer
    xlabel('x_t');
    ylabel('y_t');
    legend('TO reference path (trailer)','actual trailer','Location','Best');
    title('Trailer path: TO reference vs actual');
end
%}
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
