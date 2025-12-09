function log = simulate_towing_sample_hold(controller, tspan, params)
% controller : @(x,t) -> [F_drive, tau_r] (两个输出)
% tspan      : 1x(N+1) or (N+1)x1

    tspan = tspan(:);
    N = numel(tspan);
    x0 = params.x0(:);

    nx = numel(x0);
    X  = zeros(N, nx);
    U  = zeros(N, 2);

    X(1,:) = x0.';
    x = x0;

    for k = 1:N-1
        t  = tspan(k);
        dt = tspan(k+1) - tspan(k);
        if any(isnan(x))
        error('NaN in state BEFORE controller at step k=%d, t=%.4f', k, t);
    end
        % 控制器：显式两个输出
        [F_drive, tau_r] = controller(x, t);
        u = [F_drive; tau_r];
        u = u(:);
        U(k,:) = u.';

        % RK4 积分一步
        f = @(xx) towing_dynamics_full(xx, u, params);

          % === 检查每一个 k1,k2,k3,k4 是否出现 NaN ===
    k1 = f(x);
    if any(isnan(k1))
        error('NaN in k1 at step k=%d, t=%.4f', k, t);
    end

    x_mid = x + 0.5*dt*k1;
    k2    = f(x_mid);
    if any(isnan(k2))
        % 先把出事时刻的东西存下来
        dbg.x      = x;
        dbg.k1     = k1;
        dbg.x_mid  = x_mid;
        dbg.u      = u;
        dbg.t      = t;
        dbg.dt     = dt;
        dbg.k_step = k;
        save('debug_nan_k2.mat','dbg');

        error('NaN in k2 at step k=%d, t=%.4f (saved debug_nan_k2.mat)', k, t);
    end

    k3 = f(x + 0.5*dt*k2);
    if any(isnan(k3))
        error('NaN in k3 at step k=%d, t=%.4f', k, t);
    end

    k4 = f(x + dt*k3);
    if any(isnan(k4))
        error('NaN in k4 at step k=%d, t=%.4f', k, t);
    end

        x = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
         if any(isnan(x))
        error('NaN in state AFTER update at step k=%d, t=%.4f', k, t);
    end
        X(k+1,:) = x.';
    end

    U(N,:) = U(N-1,:);

    log.t = tspan;
    log.X = X;
    log.u = U;
end
