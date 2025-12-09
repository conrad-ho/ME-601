function [x_sol, u_sol, to_dbg] = towing_trajopt_hitch_lambda(dt, N, ref_to, x0, params)
% towing_trajopt_hitch_lambda
%   multiple shooting TO for robot+trailer.
%   目标：只对 hitch 进行轨迹优化，同时在 cost 中考虑
%        - trailer 对 hitch 的约束反作用力 lambda
%        - 铰接角 phi 的大小
%        - robot 朝向相对路径切向的偏差（yaw tracking）
%
% inputs:
%   dt      – time step
%   N       – horizon length (number of control steps, states are N+1)
%   ref_to  – struct with reference info:
%             .t    : 1 x (N+1) time grid (optional)
%             .p_h  : 2 x (N+1) hitch reference positions [xh_ref; yh_ref]
%                     (兼容 .p_hitch / .p_tr 作为 fallback)
%   x0      – initial state (nx x 1)
%   params  – struct with model parameters (masses, lengths, etc.)
%
% outputs:
%   x_sol   – (N+1) x nx optimal state trajectory
%   u_sol   – N x nu optimal control trajectory
%   to_dbg  – struct with debug info

    import casadi.*

    %% ============ 解析 hitch 参考轨迹 ============ 
    if isfield(ref_to, 'p_h')
        p_hitch_ref = ref_to.p_h;          % 推荐：显式给 hitch path
    elseif isfield(ref_to, 'p_hitch')
        p_hitch_ref = ref_to.p_hitch;      % 兼容旧字段名
    elseif isfield(ref_to, 'p_tr')
        warning('towing_trajopt_hitch_lambda: using ref_to.p_tr as hitch reference (no p_h/p_hitch).');
        p_hitch_ref = ref_to.p_tr;
    else
        error('ref_to must contain field p_h (or p_hitch / p_tr) of size 2 x (N+1).');
    end

    % === 从 hitch path 预计算路径切向 yaw_ref（running + terminal cost 都要用） ===
    yaw_ref = zeros(1, N+1);
    for k = 1:N
        dp = p_hitch_ref(:,k+1) - p_hitch_ref(:,k);   % 2x1
        yaw_ref(k) = atan2(dp(2), dp(1));
    end
    yaw_ref(N+1) = yaw_ref(N);

    %% dimensions and bounds
    nx = 12;    % [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
    nu = 2;     % [F_drive; tau_r]

    % control bounds
    if isfield(params, 'Fmax')
        Fmax = params.Fmax;
    else
        Fmax = 1e3;
    end
    if isfield(params, 'Taumax')
        Taumax = params.Taumax;
    else
        Taumax = 1e3;
    end
    u_min = [-Fmax; -Taumax];
    u_max = [ Fmax;  Taumax];

    % state bounds（先放宽）
    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% symbolic vars
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    % dynamics with lambda
    [xdot_sym, lambda_sym] = towing_dynamics_full_with_lambda(x_sym, u_sym, params);
    f_dyn = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});
    f_lam = Function('f_lam', {x_sym, u_sym}, {lambda_sym});

    %% multiple shooting vars
    X = SX.sym('X', nx, N+1);
    U = SX.sym('U', nu, N);

    %% cost 权重（优先级：hitch > yaw / phi / lambda > control）
    % 1) hitch tracking（主任务）
    Qh   = diag([30, 30]);       % running
    Qh_f = diag([100, 100]);     % terminal

    % 2) 铰接角 phi
    w_phi   = 50;                % running
    w_phi_f = 200;               % terminal

    % 3) 约束力 lambda（trailer 对 hitch 的反作用力）
    % lambda = [lambda_holo_x; lambda_holo_y; lambda_nonholo]
    w_lam_xy      = 1e-3;        % 两个 holonomic 分量权重
    w_lam_nonholo = 1e-3;        % nonholonomic 分量权重
    W_lambda      = diag([w_lam_xy, w_lam_xy, w_lam_nonholo]);

    % 4) 控制正则
    R  = diag([0.1, 1.0]);       % u magnitude
    S  = diag([0.05, 0.05]);     % u smoothness

    % 5) yaw 正则：yaw tracking + yaw rate
    w_wr   = 1.5;               % yaw rate 正则 (wr^2)
    w_yaw  = 100;                % running yaw tracking (wrap(θr-ψref))^2
    w_yaw_f = 5.0;               % 终端 yaw 对齐路径方向

    %% objective & constraints
    obj = 0;
    g   = [];

    % 初始状态约束
    g = [g; X(:,1) - x0(:)];

    %% main loop
    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        xkp = X(:,k+1);

        % --- dynamics (Euler) ---
        fk  = f_dyn(xk, uk);
        xk_euler = xk + dt * fk;
        g = [g; xkp - xk_euler];

        % --- 1) hitch tracking ---
        p_hitch_ref_k = p_hitch_ref(:,k);
        p_hitch_k     = hitch_from_state(xk, params);
        e_hitch_k     = p_hitch_k - p_hitch_ref_k;
        obj           = obj + e_hitch_k.' * Qh * e_hitch_k;

        % --- 2) yaw tracking: wrap(θ_r - ψ_ref) 到 [-pi, pi] ---
        theta_r_k = xk(3);
        psi_ref_k = yaw_ref(k);  % scalar double, 当常数
        e_yaw_k   = atan2( sin(theta_r_k - psi_ref_k), ...
                           cos(theta_r_k - psi_ref_k) );
        obj       = obj + w_yaw * (e_yaw_k^2);

        % --- 3) yaw rate 正则 ---
        wr_k = xk(6);
        obj  = obj + w_wr * (wr_k^2);

        % --- 4) 铰接角 phi penalty ---
        phi_k = hitch_angle_from_state(xk);
        obj   = obj + w_phi * (phi_k^2);

        % --- 5) lambda penalty（trailer 反作用力） ---
        lambda_k = f_lam(xk, uk);          % 3x1
        obj      = obj + lambda_k.' * W_lambda * lambda_k;

        % --- 6) control magnitude ---
        obj = obj + uk.' * R * uk;

        % --- 7) control smoothness ---
        if k > 1
            ukm1 = U(:,k-1);
            duk  = uk - ukm1;
            obj  = obj + duk.' * S * duk;
        end
    end

    %% terminal cost
    xN = X(:,N+1);

    % hitch terminal tracking
    p_hitch_ref_N = p_hitch_ref(:,N+1);
    p_hitch_N     = hitch_from_state(xN, params);
    e_hitch_N     = p_hitch_N - p_hitch_ref_N;
    obj           = obj + e_hitch_N.' * Qh_f * e_hitch_N;

    % phi terminal penalty
    phi_N = hitch_angle_from_state(xN);
    obj   = obj + w_phi_f * (phi_N^2);

    % 终端 yaw 对齐路径方向 (也做 wrap)
    theta_r_N = xN(3);
    psi_ref_N = yaw_ref(N+1);
    e_yaw_N   = atan2( sin(theta_r_N - psi_ref_N), ...
                       cos(theta_r_N - psi_ref_N) );
    obj       = obj + w_yaw_f * (e_yaw_N^2);

    %% pack decision variables
    OPT_vars = [reshape(X, nx*(N+1), 1);
                reshape(U, nu*N, 1)];

    %% build NLP
    nlp = struct;
    nlp.x = OPT_vars;
    nlp.f = obj;
    nlp.g = g;

    opts = struct;
    opts.ipopt.print_level = 3;
    opts.print_time = 0;
    opts.ipopt.acceptable_tol = 1e-6;
    opts.ipopt.acceptable_obj_change_tol = 1e-4;

    solver = nlpsol('solver', 'ipopt', nlp, opts);

    %% bounds on constraints (全部是等式，0)
    ng  = length(g);
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

    %% bounds on variables
    nX = nx*(N+1);
    nU = nu*N;

    lbx = -inf(size(OPT_vars));
    ubx =  inf(size(OPT_vars));

    lbx(1:nX) = repmat(x_min, N+1, 1);
    ubx(1:nX) = repmat(x_max, N+1, 1);

    lbx(nX+1:nX+nU) = repmat(u_min, N, 1);
    ubx(nX+1:nX+nU) = repmat(u_max, N, 1);

    %% initial guess
    X_init = repmat(x0(:), 1, N+1);
    U_init = zeros(nu, N);

    x0_guess = [reshape(X_init, nX, 1);
                reshape(U_init, nU, 1)];

    %% solve
    sol = solver('x0', x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    % unpack
    X_opt = reshape(sol_vec(1:nX), nx, N+1);
    U_opt = reshape(sol_vec(nX+1:end), nu, N);

    x_sol = X_opt.';   % (N+1) x nx
    u_sol = U_opt.';   % N x nu

    % debug info
    to_dbg = struct;
    to_dbg.obj      = full(sol.f);
    to_dbg.sol_vec  = sol_vec;
    to_dbg.stats    = solver.stats();

end

%% ======== helper functions (保持你现有的定义) ========

function xdot = towing_dynamics_full(x, u, params)
    import casadi.*

    q = [x(1:3); x(7:9)];
    v = [x(4:6); x(10:12)];

    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;
    B    = dyn.B;
    J    = dyn.J;
    dotJ = dyn.dotJ;
    vgen = dyn.v;

    K   = [M, -J';
           J, SX.zeros(size(J,1), size(J,1))];
    rhs = [B*u;
          -dotJ*vgen];

    sol = K \ rhs;
    a   = sol(1:6);

    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end

function p_h = hitch_from_state(x,params)
    d  = params.d;
    xr = x(1);  yr = x(2);  thetar = x(3);
    p_h = [xr - d*cos(thetar);
           yr - d*sin(thetar)];
end

function phi = hitch_angle_from_state(x)
    theta_r = x(3);
    theta_t = x(9);
    phi = theta_t - theta_r;
end
