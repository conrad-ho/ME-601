function [x_sol, u_sol, to_dbg] = towing_trajopt(dt, N, ref_to, x0, params)
% towing_trajopt
%   basic multiple shooting trajectory optimization for a robot+trailer system.
%   trailer is the main tracking target, robot has weaker soft constraints.
%
% inputs:
%   dt      – time step
%   N       – horizon length (number of control steps, states are N+1)
%   ref_to 结构：
%   ref_to.t        : 1 x (N+1) 时间网格（可选）
%   ref_to.p_hitch  : 2 x (N+1) hitch 期望位置 [xh_ref; yh_ref]
%   x0      – initial state (nx x 1)
%   params  – struct with model parameters (masses, lengths, etc.)
%
% outputs:
%   x_sol   – (N+1) x nx optimal state trajectory
%   u_sol   – N x nu optimal control trajectory
%   to_dbg  – struct with debug info (solver stats, cost terms, etc.)

    import casadi.*
    % ============ parse reference path (hitch path) ============
    % 推荐：ref.p_hitch，兼容旧写法：ref.p_tr
    if isfield(ref_to, 'p_hitch')
        p_hitch_ref = ref_to.p_hitch;      % 2 x (N+1)
    elseif isfield(ref_to, 'p_tr')
        p_hitch_ref = ref_to.p_tr;         % 2 x (N+1)，旧字段名
    else
        error('ref must contain field p_hitch or p_tr (2 x (N+1)).');
    end

    %% dimensions and basic bounds
    % TODO: adjust nx, nu 根据你实际的状态和控制维度来改
    nx = 12;    % [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
    nu = 2;   % example: [F, tau]

    % control bounds (example)
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
    u_min = [-Fmax; -Taumax];  % 2x1
    u_max = [ Fmax;  Taumax];  % 2x1

    % optional state bounds (可以先不严，后面再收紧)
    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% symbolic variables: state, control
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    % ========= dynamics =========
    xdot_sym = towing_dynamics_full(x_sym, u_sym, params);

    f_dyn = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});

    %% multiple shooting variables
    X = SX.sym('X', nx, N+1);   % states over horizon
    U = SX.sym('U', nu, N);     % controls over horizon

    %% cost weights
    % trailer path tracking
    Qp_tr = diag([10, 10]);     % trailer position tracking weight

    % robot soft reference
    Qp_r  = diag([1, 1]);       % robot soft position weight
    w_r   = 0.1;                % scalar factor for robot soft penalty

    % trailer hitch angle penalty
    w_phi    = 1.0;             % running penalty
    w_phi_f  = 5.0;             % terminal penalty

    % control effort and smoothness
    R  = diag([0.1, 0.1]);      % control magnitude weight
    S  = diag([1.0, 1.0]);      % control difference weight (smoothness)

    % terminal trailer tracking weight
    Qp_tr_f = diag([50, 50]);

    %% objective and constraints
    obj = 0;
    g   = [];

    % initial condition constraint: X(:,1) = x0
    g = [g; X(:,1) - x0(:)];

    %% main loop over horizon
    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        xkp = X(:,k+1);  % x_{k+1}

        % ---- dynamics constraint: x_{k+1} = x_k + dt * f(x_k, u_k) ----
        fk  = f_dyn(xk, uk);
        xk_euler = xk + dt * fk;
        g = [g; xkp - xk_euler];

        % ---- trailer position tracking ----
        % trailer ref at step k
        p_tr_ref_k = ref_to.p_tr(:,k);      % 2x1 TO's reference path
        p_tr_k     =  hitch_from_state(xk, params); % actual x_r read from xk

        e_tr = p_tr_k - p_tr_ref_k;
        obj  = obj + e_tr.' * Qp_tr * e_tr;
        
        % ---- robot soft reference (optional) ----
        %{
        %TODO after make TO work
        if isfield(ref, 'p_r')
            p_r_ref_k = ref.p_r(:,k);    % 2x1
        else
            % 如果没有单独给 robot 参考，可以取近似：
            % 例如和 trailer ref 一样，或者加一点几何偏移
            p_r_ref_k = p_tr_ref_k;
        end
        
        p_r_k = robot_pos_from_state(xk);   % TODO: 根据你的状态定义修改
        e_r   = p_r_k - p_r_ref_k;
        obj   = obj + w_r * (e_r.' * Qp_r * e_r);
        %}

        % ---- trailer hitch angle penalty ---- 等work以后调参
        %{
        phi_k = hitch_angle_from_state(xk); % TODO: 根据你的状态定义修改
        obj   = obj + w_phi * (phi_k^2);
        %}
        % ---- control magnitude ----
        obj   = obj + uk.' * R * uk;

        % ---- control smoothness penalty: ||u_k - u_{k-1}||_S^2 ----
        %{
        %等work以后调整
        if k > 1
            ukm1 = U(:,k-1);
            duk  = uk - ukm1;
            obj  = obj + duk.' * S * duk;
        end
        %}
    end

    %% terminal cost at k = N+1
    xN = X(:,N+1);

    p_tr_ref_N = ref_to.p_tr(:,N+1);
    p_tr_N     = hitch_from_state(xN,params);
    e_tr_N     = p_tr_N - p_tr_ref_N;
    obj        = obj + e_tr_N.' * Qp_tr_f * e_tr_N;
    
    phi_N      = hitch_angle_from_state(xN);
    obj        = obj + w_phi_f * (phi_N^2);

    %% pack decision variables
    % decision vector: [vec(X); vec(U)]
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
    % equality constraints: initial condition + dynamics
    % g has (nx*(N+1)) entries from:
    %   - nx from X(:,1)-x0
    %   - nx*N from dynamics
    % 没有别的不等式约束的话，全部是等式 => lbg=ubg=0
    ng = nx*(N+1);
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

    %% bounds on variables
    nX = nx*(N+1);
    nU = nu*N;

    lbx = -inf(size(OPT_vars));
    ubx =  inf(size(OPT_vars));

    % state bounds
    lbx(1:nX) = repmat(x_min, N+1, 1);
    ubx(1:nX) = repmat(x_max, N+1, 1);

    % control bounds
    lbx(nX+1:nX+nU) = repmat(u_min, N, 1);
    ubx(nX+1:nX+nU) = repmat(u_max, N, 1);

    %% initial guess
    % TODO: 这里可以用 path interpolation / simple rollout 给一个更聪明的初值
    X_init = repmat(x0(:), 1, N+1);  % 先所有步都用初始状态
    U_init = zeros(nu, N);           % 控制都从零开始

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

%% ======== helper functions (placeholders) ========

function xdot = towing_dynamics_full(x, u, params)
    import casadi.*

    % 拆 q,v
    q = [x(1:3); x(7:9)];
    v = [x(4:6); x(10:12)];

    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;      % 6x6
    B    = dyn.B;      % 6x2
    J    = dyn.J;      % 3x6
    dotJ = dyn.dotJ;   % 3x6
    vgen = dyn.v;      % 6x1

    % KKT 系统
    K   = [M, -J';
           J, SX.zeros(size(J,1), size(J,1))];    % 9x9
    rhs = [B*u;
          -dotJ*vgen];                            % 9x1

    sol = K \ rhs;
    a   = sol(1:6);   % generalized acceleration

    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end


function p_tr = hitch_from_state(x,params)
    d  = params.d;
    xr = x(1);  yr = x(2);  thetar = x(3);

    % hitch 点在世界坐标下的位置
    p_tr = [xr - d*cos(thetar);
            yr - d*sin(thetar)];
end

function p_r = robot_pos_from_state(x)
    % TODO: modify according to your state definition
    % example: x = [xr; yr; thetar; xt; yt; thetat]
    xr = x(1);
    yr = x(2);
    p_r = [xr; yr];     
end

function phi = hitch_angle_from_state(x)
    % TODO: modify according to your hitch definition
    % example: relative yaw between trailer and robot
    theta_r = x(3);
    theta_t = x(9);
    phi = theta_t - theta_r;
end
