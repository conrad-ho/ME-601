function debug_yaw_sim(log_file)
% debug_yaw_sim('log.mat')

    data   = load(log_file);
    log    = data.log;
    ref_to = data.ref_to;

    t_log = log.t;           % N_log x 1
    X     = log.X;           % N_log x nx

    % ---- hitch path & yaw_ref on TO grid ----
    p_h = ref_to.p_hitch;         % 2 x (N_to+1)
    Nto = size(p_h, 2) - 1;
    t_to = linspace(t_log(1), t_log(end), Nto+1);  % 这里按你 TO dt 改也行

    yaw_ref_to = zeros(1, Nto+1);
    for k = 1:Nto
        dp = p_h(:,k+1) - p_h(:,k);
        yaw_ref_to(k) = atan2(dp(2), dp(1));
    end
    yaw_ref_to(Nto+1) = yaw_ref_to(Nto);

    % 插值到仿真时间上
    yaw_ref_log = interp1(t_to, yaw_ref_to, t_log, 'linear', 'extrap');

    % ---- sim 中的 robot yaw ----
    theta_r_log = X(:,3);    % 按你的状态定义改
    theta_r_unwrap = unwrap(theta_r_log);

    % ---- yaw 误差 ----
    e_raw  = theta_r_log - yaw_ref_log;
    e_wrap = atan2( sin(e_raw), cos(e_raw) );

    figure;
    subplot(3,1,1);
    plot(t_log, yaw_ref_log, '-', t_log, theta_r_log, '--');
    legend('\psi_{ref}(t)','\theta_r(t)');
    grid on; title('Yaw in simulation (raw)');

    subplot(3,1,2);
    plot(t_log, e_raw); grid on;
    ylabel('e\_raw');

    subplot(3,1,3);
    plot(t_log, e_wrap); grid on;
    ylabel('e\_wrap');
    xlabel('time');
    title('e\_wrap = atan2(sin(..),cos(..))');
end
