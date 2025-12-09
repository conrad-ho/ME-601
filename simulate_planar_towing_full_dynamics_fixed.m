function log = simulate_planar_towing_full_dynamics_fixed(controller_handle, tspan, params)
%SIMULATE_PLANAR_TOWING_FULL_DYNAMICS_FIXED
%   Fixed-step RK4 simulation for robot+trailer full dynamics.
%
%   controller_handle : @(x,t) -> [F_drive, tau_r, dbg]
%   tspan             : [t0, tf]
%   params.x0         : 12x1 initial state
%
%   log.t  [M x 1]
%   log.X  [M x 12]
%   log.u  [M x 2]
%
%   如果状态范数超过某个上限（比如 1e3），仿真会提前停止，避免 NaN。

    if ~isfield(params,'x0')
        error('simulate_fixed:missing_x0', 'params.x0 (12x1) is required.');
    end

    x0 = params.x0(:);
    if numel(x0) ~= 12
        error('simulate_fixed:bad_x0_dim', ...
              'params.x0 must be 12x1, got size %s', mat2str(size(x0)));
    end

    t0 = tspan(1);
    tf = tspan(end);

    if isfield(params,'dt_sim')
        dt = params.dt_sim;
    else
        dt = 1e-3;
    end

    N = floor((tf - t0)/dt);
    t_vec = t0 + (0:N)'*dt;

    X = zeros(N+1, 12);
    U = zeros(N+1, 2);

    X(1,:) = x0.';

    % 状态幅度上限（防止跑到 1e40 那种鬼畜位置）
    X_max_norm = 1e3;

    M_used = N+1;  % 实际用到的步数

    for k = 1:N
        t  = t_vec(k);
        xk = X(k,:).';

        % 若状态已经太大，提前终止
        if norm(xk) > X_max_norm
            fprintf(2, 'simulate_fixed: state norm %.3e > %.1e at step %d, t=%.3f, stop.\n', ...
                    norm(xk), X_max_norm, k, t);
            M_used = k;
            break;
        end

        % controller
        [F_drive, tau_r, ~] = controller_handle(xk, t);
        u = [F_drive; tau_r];
        U(k,:) = u.';

        f = @(xx) towing_dynamics_full(xx, u, params);

        try
            k1 = f(xk);
            k2 = f(xk + 0.5*dt*k1);
            k3 = f(xk + 0.5*dt*k2);
            k4 = f(xk + dt*k3);
        catch ME
            fprintf(2, 'Dynamics crashed at step k = %d, t = %.6f\n', k, t);
            fprintf(2, '  message: %s\n', ME.message);
            fprintf(2, '  xk'' = %s\n', mat2str(xk(:).', 6));
            fprintf(2, '  u''  = %s\n', mat2str(u(:).', 6));
            rethrow(ME);
        end

        x_next = xk + dt/6 * (k1 + 2*k2 + 2*k3 + k4);

        if any(~isfinite(x_next))
            error('simulate_fixed:NaN_in_state', ...
                  'NaN/Inf in state at step %d, t = %.6f', k, t);
        end

        X(k+1,:) = x_next.';
    end

    % 截断到实际用到的长度
    X = X(1:M_used, :);
    U = U(1:M_used, :);
    t_vec = t_vec(1:M_used);

    log = struct();
    log.t = t_vec;
    log.X = X;
    log.u = U;
end
