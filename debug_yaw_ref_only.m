function debug_yaw_ref_only(ref_file)
% debug_yaw_ref_only('ref_to.mat')  % or load from log.mat

    data = load(ref_file);   % 里面要有 ref_to
    ref_to = data.ref_to;

    p_h = ref_to.p_hitch;    % 2 x (N+1)  [xh; yh]
    N   = size(p_h, 2) - 1;

    yaw_ref = zeros(1, N+1);
    for k = 1:N
        dp = p_h(:,k+1) - p_h(:,k);
        yaw_ref(k) = atan2(dp(2), dp(1));
    end
    yaw_ref(N+1) = yaw_ref(N);

    % unwrap 版本（看是不是有 2pi 的跳）
    yaw_ref_unwrap = unwrap(yaw_ref);

    figure;
    subplot(3,1,1);
    plot(p_h(1,:), p_h(2,:), 'o-');
    axis equal; grid on;
    title('Hitch reference path');

    subplot(3,1,2);
    plot(0:N, yaw_ref, 'o-'); grid on;
    ylabel('\psi_{ref} (raw)');
    title('raw yaw\_ref (atan2)');

    subplot(3,1,3);
    plot(0:N, yaw_ref_unwrap, 'o-'); grid on;
    ylabel('\psi_{ref} unwrapped');
    title('unwrap(yaw\_ref)');
end
