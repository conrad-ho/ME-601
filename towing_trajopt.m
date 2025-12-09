function [x_sol, u_sol, to_dbg] = towing_trajopt(dt, N, ref_to, x0, params)
% towing_trajopt
%   multiple shooting TO for robot+trailer.
%   优先级：
%     1) hitch 强 tracking
%     2) 铰接角 phi 强约束（避免 trailer 暴走）
%     3) trailer COM 中等权重 tracking
%
% inputs:
%   dt      – time step
%   N       – horizon length (number of control steps, states are N+1)
%   ref_to  – struct with reference info:
%             .t     : 1 x (N+1) time grid (optional)
%             .p_h   : 2 x (N+1) hitch reference positions [xh_ref; yh_ref]
%             .p_tr  : 2 x (N+1) trailer COM reference [xt_ref; yt_ref] (可选)
%   x0      – initial state (nx x 1)
%   params  – struct with model parameters (masses, lengths, etc.)
%
% outputs:
%   x_sol   – (N+1) x nx optimal state trajectory
%   u_sol   – N x nu optimal control trajectory
%   to_dbg  – struct with debug info

    import casadi.*
    % ============ parse hitch reference path ============
    if isfield(ref_to, 'p_h')
        p_hitch_ref = ref_to.p_h;          % 2 x (N+1) hitch path
    elseif isfield(ref_to, 'p_hitch')
        p_hitch_ref = ref_to.p_hitch;      % 兼容旧字段名
    elseif isfield(ref_to, 'p_tr')
        warning('towing_trajopt: using ref_to.p_tr as hitch reference (no p_h/p_hitch field).');
        p_hitch_ref = ref_to.p_tr;
    else
        error('ref_to must contain field p_h (or p_hitch / p_tr) of size 2 x (N+1).');
    end

    % === 从 hitch path 预计算期望行驶方向 yaw_ref ===
    yaw_ref = zeros(1, N+1);
    for k = 1:N
        dp = p_hitch_ref(:,k+1) - p_hitch_ref(:,k);   % 2x1
        yaw_ref(k) = atan2(dp(2), dp(1));             % 路径切线方向
    end
    yaw_ref(N+1) = yaw_ref(N);  % 末端沿用最后一个方向

    %% dimensions and basic bounds
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

    % state bounds
    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% symbolic variables: state, control
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    % ========= dynamics =========
    xdot_sym = towing_dynamics_full(x_sym, u_sym, params);
    f_dyn = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});

    %% multiple shooting variables
    X = SX.sym('X', nx, N+1);   % states
    U = SX.sym('U', nu, N);     % controls

    %% cost weights (三层优先级)
    % 1) hitch 强 tracking
    Qh    = diag([30, 30]);      % running
    Qh_f  = diag([100, 100]);    % terminal

    % 2) 铰接角 phi 强约束
    w_phi   = 50;                 % running
    w_phi_f = 200;                % terminal

    % 3) trailer COM 中等权重 (相对 hitch 小一些)
    % 只有在 ref_to.p_tr 存在时才启用
    Qtr    = 0.3 * Qh;            % ~ 0.3 * hitch 权重
    Qtr_f  = 0.4 * Qh_f;

    % control 正则
    R  = diag([0.1, 0.1]);        % control magnitude
    S  = diag([0.05, 0.05]);      % control smoothness

    % yaw 相关
    w_wr    = 0.01;               % yaw rate 正则
    w_yaw_f = 0.10;               % 终端 yaw 对齐路径方向

    %% objective and constraints
    obj = 0;
    g   = [];

    % initial condition
    g = [g; X(:,1) - x0(:)];

    %% main loop
    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        xkp = X(:,k+1);

        % ------ dynamics (Euler) ------
        fk = f_dyn(xk, uk);
        xk_euler = xk + dt * fk;
        g = [g; xkp - xk_euler];

        % ------ 1) hitch tracking (主目标) ------
        p_hitch_ref_k = p_hitch_ref(:,k);
        p_hitch_k     = hitch_from_state(xk, params);
        e_hitch_k     = p_hitch_k - p_hitch_ref_k;
        obj           = obj + e_hitch_k.' * Qh * e_hitch_k;

        % ------ yaw-rate 正则 ------
        wr_k = xk(6);                 % robot yaw rate
        obj  = obj + w_wr * (wr_k^2);

        % ------ 3) trailer COM 软 tracking (中等优先级) ------
        if isfield(ref_to, 'p_tr')
            p_tr_com_k     = [xk(7); xk(8)];       % trailer COM
            p_tr_com_ref_k = ref_to.p_tr(:,k);     % reference
            e_tr_com_k     = p_tr_com_k - p_tr_com_ref_k;
            obj            = obj + e_tr_com_k.' * Qtr * e_tr_com_k;
        end

        % ------ 2) 铰接角 phi penalty (避免暴走) ------
        phi_k = hitch_angle_from_state(xk);
        obj   = obj + w_phi * (phi_k^2);

        % ------ control magnitude ------
        obj = obj + uk.' * R * uk;

        % ------ control smoothness ------
        if k > 1
            ukm1 = U(:,k-1);
            duk  = uk - ukm1;
            obj  = obj + duk.' * S * duk;
        end
    end

    %% terminal cost
    xN = X(:,N+1);

    % hitch 终端 tracking
    p_hitch_ref_N = p_hitch_ref(:,N+1);
    p_hitch_N     = hitch_from_state(xN, params);
    e_hitch_N     = p_hitch_N - p_hitch_ref_N;
    obj           = obj + e_hitch_N.' * Qh_f * e_hitch_N;

    % 铰接角终端 penalty
    phi_N = hitch_angle_from_state(xN);
    obj   = obj + w_phi_f * (phi_N^2);

    % trailer COM 终端 tracking
    if isfield(ref_to, 'p_tr')
        p_tr_com_N     = [xN(7); xN(8)];
        p_tr_com_ref_N = ref_to.p_tr(:,N+1);
        e_tr_com_N     = p_tr_com_N - p_tr_com_ref_N;
        obj            = obj + e_tr_com_N.' * Qtr_f * e_tr_com_N;
    end

    % 终端 yaw 对齐路径方向
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

    %% bounds on constraints
    ng  = nx*(N+1);
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

    %% solve NLP
    sol = solver('x0', x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    % unpack solution
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
