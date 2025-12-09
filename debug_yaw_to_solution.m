function debug_yaw_to_solution(to_file,log_file)
% debug_yaw_to_solution('TO_Output.mat')

    data   = load(to_file);
    logdata= load(log_file);
    x_to   = data.x_to;     % (N+1) x nx
    ref_to = logdata.ref_to;
    
    % ---- hitch path → yaw_ref（和 TO 里代码保持一致）----
    p_h = ref_to.p_hitch;    % 2 x (N+1)
    N   = size(p_h, 2) - 1;

    yaw_ref = zeros(1, N+1);
    for k = 1:N
        dp = p_h(:,k+1) - p_h(:,k);
        yaw_ref(k) = atan2(dp(2), dp(1));
    end
    yaw_ref(N+1) = yaw_ref(N);
    yaw_ref_unwrap = unwrap(yaw_ref);

    % ---- robot yaw from TO solution ----
    theta_r = x_to(:,3)';    % 按你的状态定义改列号
    theta_r_unwrap = unwrap(theta_r);

    % ---- yaw 误差：raw vs wrap ----
    e_raw  = theta_r - yaw_ref;
    e_wrap = atan2( sin(theta_r - yaw_ref), ...
                    cos(theta_r - yaw_ref) );

    figure;
    subplot(4,1,1);
    plot(p_h(1,:), p_h(2,:), 'o-'); axis equal; grid on;
    title('Hitch reference path');

    subplot(4,1,2);
    plot(0:N, yaw_ref_unwrap, 'o-', 0:N, theta_r_unwrap, 'x-');
    legend('\psi_{ref} unwrap','\theta_r unwrap');
    grid on; title('Yaw (ref vs TO solution)');

    subplot(4,1,3);
    plot(0:N, e_raw, 'o-'); grid on;
    ylabel('e\_yaw raw');
    title('e\_raw = \theta_r - \psi_{ref}');

    subplot(4,1,4);
    plot(0:N, e_wrap, 'o-'); grid on;
    ylabel('e\_wrap');
    title('e\_wrap = atan2(sin(..),cos(..))');
end
