function compare_to_tuning(file_before, file_after, params)
% compare_to_tuning
%   load two logs (before / after tuning to) and compare performance.
%
% inputs:
%   file_before : string, e.g. 'log_before.mat'
%   file_after  : string, e.g. 'log_after.mat'
%   params.d    : hitch offset (robot com to hitch)

    if ~isfield(params,'d')
        error('params.d (hitch offset) is required');
    end

    % load both runs
    [log_b, ref_b] = load_one_run(file_before);
    [log_a, ref_a] = load_one_run(file_after);

    % compute metrics for both
    metrics_b = compute_metrics(log_b, ref_b, params);
    metrics_a = compute_metrics(log_a, ref_a, params);

    % print summary to command window
    print_metrics_table(metrics_b, metrics_a);

    % make comparison plots
    make_plots(log_b, ref_b, metrics_b, log_a, ref_a, metrics_a);
end

% -------------------------------------------------------------------------
function [log, ref_to] = load_one_run(fname)
    data = load(fname);

    if ~isfield(data,'log')
        error('%s does not contain variable ''log''.', fname);
    end
    if ~isfield(data,'ref_to')
        error('%s does not contain variable ''ref_to''.', fname);
    end

    log    = data.log;
    ref_to = data.ref_to;
end

% -------------------------------------------------------------------------
function m = compute_metrics(log, ref_to, params)
% compute all numeric metrics for one run

    % basic alignment
    t = log.t(:);
    N_log = length(t);
    p_tr_ref = ref_to.p_tr.';    % 2xN -> Nx2
    p_ht_ref = ref_to.p_hitch.'; % 2xN -> Nx2

    N_ref = size(p_tr_ref,1);
    N = min(N_log, N_ref);

    t = t(1:N);
    y_tr = log.y(1:N,:);        % trailer com actual (x_t, y_t)
    p_tr = p_tr_ref(1:N,:);     % trailer com ref
    p_ht = p_ht_ref(1:N,:);     % hitch ref

    % trailer tracking error
    e_tr = y_tr - p_tr;
    e_tr_norm = vecnorm(e_tr,2,2);
    m.trailer_err_max = max(e_tr_norm);
    m.trailer_err_rms = sqrt(mean(e_tr_norm.^2));

    % if X available, compute hitch actual and robot/trailer yaw info
    if isfield(log,'X') && size(log.X,2) >= 9
        X = log.X(1:N,:);
        xr = X(:,1);
        yr = X(:,2);
        thetar = X(:,3);
        xt = X(:,7);
        yt = X(:,8);
        thetat = X(:,9);

        d = params.d;

        % hitch actual position from robot com
        xh = xr - d.*cos(thetar);
        yh = yr - d.*sin(thetar);
        y_hitch = [xh, yh];

        % hitch tracking error
        e_ht = y_hitch - p_ht;
        e_ht_norm = vecnorm(e_ht,2,2);
        m.hitch_err_max = max(e_ht_norm);
        m.hitch_err_rms = sqrt(mean(e_ht_norm.^2));

        % yaw behaviour
        thetar_un = unwrap(thetar);
        thetat_un = unwrap(thetat);

        m.robot_yaw_net_change    = thetar_un(end) - thetar_un(1);
        m.robot_yaw_total_swing   = max(thetar_un) - min(thetar_un);
        m.robot_yaw_std           = std(thetar_un);

        m.trailer_yaw_net_change  = thetat_un(end) - thetat_un(1);
        m.trailer_yaw_total_swing = max(thetat_un) - min(thetat_un);
        m.trailer_yaw_std         = std(thetat_un);
    else
        m.hitch_err_max  = NaN;
        m.hitch_err_rms  = NaN;
        m.robot_yaw_net_change    = NaN;
        m.robot_yaw_total_swing   = NaN;
        m.robot_yaw_std           = NaN;
        m.trailer_yaw_net_change  = NaN;
        m.trailer_yaw_total_swing = NaN;
        m.trailer_yaw_std         = NaN;
    end

    % control effort / smoothness
    if isfield(log,'u')
        u  = log.u(1:N,:);
        dt = mean(diff(t)); % assume roughly uniform

        m.u_L2 = sum(sum(u.^2,2))*dt;              % integral of ||u||^2
        m.F_max = max(abs(u(:,1)));
        m.tau_max = max(abs(u(:,2)));

        du = diff(u);
        m.du_rms = sqrt(mean(sum(du.^2,2)));       % rough smoothness indicator
    else
        m.u_L2   = NaN;
        m.F_max  = NaN;
        m.tau_max = NaN;
        m.du_rms = NaN;
    end

    % you可以在这里继续加入 path length, overshoot 等指标
end

% -------------------------------------------------------------------------
function print_metrics_table(mb, ma)
% print side-by-side comparison in command window

    fprintf('\n============== TO tuning comparison ==============\n');
    fprintf('                 before         after\n');
    fprintf('--------------------------------------------------\n');
    fprintf('trailer err max   %8.4f      %8.4f\n', mb.trailer_err_max, ma.trailer_err_max);
    fprintf('trailer err rms   %8.4f      %8.4f\n', mb.trailer_err_rms, ma.trailer_err_rms);
    fprintf('hitch err max     %8.4f      %8.4f\n', mb.hitch_err_max,   ma.hitch_err_max);
    fprintf('hitch err rms     %8.4f      %8.4f\n', mb.hitch_err_rms,   ma.hitch_err_rms);
    fprintf('u L2 cost         %8.4f      %8.4f\n', mb.u_L2,            ma.u_L2);
    fprintf('max |F|           %8.4f      %8.4f\n', mb.F_max,           ma.F_max);
    fprintf('max |tau|         %8.4f      %8.4f\n', mb.tau_max,         ma.tau_max);
    fprintf('du rms            %8.4f      %8.4f\n', mb.du_rms,          ma.du_rms);
    fprintf('robot yaw netΔ    %8.4f      %8.4f\n', mb.robot_yaw_net_change,    ma.robot_yaw_net_change);
    fprintf('robot yaw swing   %8.4f      %8.4f\n', mb.robot_yaw_total_swing,   ma.robot_yaw_total_swing);
    fprintf('robot yaw std     %8.4f      %8.4f\n', mb.robot_yaw_std,           ma.robot_yaw_std);
    fprintf('==================================================\n\n');
end

% -------------------------------------------------------------------------
function make_plots(log_b, ref_b, mb, log_a, ref_a, ma)
    % trailer com path comparison
    figure('Color','w'); hold on; grid on; axis equal;
    p_tr_b = ref_b.p_tr.';   % Nx2
    p_tr_a = ref_a.p_tr.';
    plot(p_tr_b(:,1), p_tr_b(:,2), 'r--', 'LineWidth', 1.0);
    plot(log_b.y(:,1), log_b.y(:,2), 'b-', 'LineWidth', 1.2);
    plot(log_a.y(:,1), log_a.y(:,2), 'g-', 'LineWidth', 1.2);
    xlabel('x_t'); ylabel('y_t');
    legend('TO ref (before)','actual trailer before',...
           'actual trailer after','Location','Best');
    title('trailer com path: before vs after tuning');

    % error norm over time (trailer)
    figure('Color','w');
    subplot(2,1,1); hold on; grid on;
    N_b = min(length(log_b.t), size(ref_b.p_tr,2));
    N_a = min(length(log_a.t), size(ref_a.p_tr,2));
    t_b = log_b.t(1:N_b);
    t_a = log_a.t(1:N_a);
    e_tr_b = log_b.y(1:N_b,:).' - ref_b.p_tr(:,1:N_b);
    e_tr_a = log_a.y(1:N_a,:).' - ref_a.p_tr(:,1:N_a);
    e_tr_b_norm = vecnorm(e_tr_b,2,1);
    e_tr_a_norm = vecnorm(e_tr_a,2,1);
    plot(t_b, e_tr_b_norm, 'b-');
    plot(t_a, e_tr_a_norm, 'g-');
    ylabel('‖e_{tr}‖');
    legend('before','after');
    title('trailer error norm');

    % control inputs comparison
    subplot(2,1,2); hold on; grid on;
    plot(log_b.t, log_b.u(:,1), 'b-', 'DisplayName','F before');
    plot(log_a.t, log_a.u(:,1), 'g-', 'DisplayName','F after');
    plot(log_b.t, log_b.u(:,2), 'b--', 'DisplayName','tau before');
    plot(log_a.t, log_a.u(:,2), 'g--', 'DisplayName','tau after');
    xlabel('t');
    ylabel('u');
    legend('Location','Best');
    title('control inputs: before vs after');

    % 你也可以另外开一个 figure 专门画 robot yaw / trailer yaw
end
