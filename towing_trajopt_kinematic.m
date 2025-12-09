function [x_sol, u_sol, to_dbg] = towing_trajopt_kinematic(dt, N, ref_to, x0_kin, params)
% towing_trajopt_kinematic
%   Reduced-DOF kinematic TO for robot+trailer:
%   state x = [xr; yr; thetar; theta_t] (4x1)
%   input u = [v_r; w_r]               (2x1)
%
%   Kinematic model:
%     dxr/dt     = v_r * cos(thetar)
%     dyr/dt     = v_r * sin(thetar)
%     dthetar/dt = w_r
%     dtheta_t/dt= (v_r / Lt) * sin(thetar - theta_t)
%
%   Geometry (robot on right, trailer on left, both initially facing +x):
%     hitch from robot COM:
%       p_h = [xr - d*cos(thetar);
%              yr - d*sin(thetar)]
%     trailer COM:
%       p_tr = p_h - [Lt*cos(theta_t);
%                     Lt*sin(theta_t)]
%
% inputs:
%   dt      : time step
%   N       : horizon length (controls N, states N+1)
%   ref_to  : struct, optional fields:
%               .t    : 1 x (N+1) time grid
%               .p_r  : 2 x (N+1) robot com reference path
%               .p_tr : 2 x (N+1) trailer com reference path
%   x0_kin  : initial kinematic state [xr0;yr0;thetar0;theta_t0]
%   params  : struct with fields
%               .d   : robot COM to hitch distance
%               .Lt  : trailer hitch to COM distance
%               .vmax, .wmax : (optional) speed bounds
%
% outputs:
%   x_sol   : (N+1) x 4 optimal kinematic state trajectory
%   u_sol   : N x 2    optimal [v_r, w_r] trajectory
%   to_dbg  : struct, debug stats

    import casadi.*

    %% dimensions
    nx = 4;   % [xr; yr; thetar; theta_t]
    nu = 2;   % [v_r; w_r]

    d  = params.d;
    Lt = params.Lt;

    %% bounds on controls
    if isfield(params,'vmax')
        vmax = params.vmax;
    else
        vmax = 2.0;   % you can tune
    end
    if isfield(params,'wmax')
        wmax = params.wmax;
    else
        wmax = 2.0;   % you can tune
    end
    u_min = [-vmax; -wmax];
    u_max = [ vmax;  wmax];

    %% state bounds (can be tightened if desired)
    x_min = [-inf; -inf; -2*pi; -2*pi];
    x_max = [ inf;  inf;  2*pi;  2*pi];

    %% symbolic variables
    x_sym = SX.sym('x', nx, 1);   % [xr; yr; thetar; theta_t]
    u_sym = SX.sym('u', nu, 1);   % [v_r; w_r]

    xr      = x_sym(1);
    yr      = x_sym(2);
    thetar  = x_sym(3);
    theta_t = x_sym(4);

    v_r = u_sym(1);
    w_r = u_sym(2);

    % kinematic towing model
    xdot_sym = [
        v_r * cos(thetar);                 % dxr/dt
        v_r * sin(thetar);                 % dyr/dt
        w_r;                               % dthetar/dt
        (v_r / Lt) * sin(thetar - theta_t) % dtheta_t/dt
    ];

    f_kin = Function('f_kin', {x_sym,u_sym}, {xdot_sym});

    %% RK4 one-step map: x_{k+1} = F_kin(x_k,u_k)
    xk_sym = SX.sym('xk', nx, 1);
    uk_sym = SX.sym('uk', nu, 1);

    k1 = f_kin(xk_sym,           uk_sym);
    k2 = f_kin(xk_sym+0.5*dt*k1, uk_sym);
    k3 = f_kin(xk_sym+0.5*dt*k2, uk_sym);
    k4 = f_kin(xk_sym+dt*k3,     uk_sym);

    xk_rk4 = xk_sym + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    F_kin  = Function('F_kin', {xk_sym,uk_sym}, {xk_rk4});

    %% multiple-shooting variables
    X = SX.sym('X', nx, N+1);   % states
    U = SX.sym('U', nu, N);     % controls

    %% build yaw_ref from robot reference path, if available
    yaw_ref = zeros(1,N+1);
    if isfield(ref_to,'p_r')
        p_r_ref = ref_to.p_r;
        for k = 1:N
            dp = p_r_ref(:,k+1) - p_r_ref(:,k);
            yaw_ref(k) = atan2(dp(2), dp(1));
        end
        yaw_ref(N+1) = yaw_ref(N);
    else
        yaw_ref(:) = 0;
    end

    %% cost weights
    % robot com tracking
    Qr    = diag([30, 30]);
    Qr_f  = diag([100, 100]);

    % trailer com tracking
    Qt    = 0.3 * Qr;
    Qt_f  = 0.4 * Qr_f;

    % articulation angle phi = theta_t - thetar
    w_phi   = 20;
    w_phi_f = 80;

    % control regularization
    R  = diag([0.1, 0.1]);    % magnitude
    S  = diag([0.05, 0.05]);  % smoothness

    % yaw alignment terminal
    w_yaw_f = 0.1;

    %% objective and constraints
    obj = 0;
    g   = [];

    % initial condition
    g = [g; X(:,1) - x0_kin(:)];

    %% main loop
    for k = 1:N
        xk = X(:,k);
        uk = U(:,k);
        xkp = X(:,k+1);

        % discrete kinematics
        xk_next = F_kin(xk,uk);
        g = [g; xkp - xk_next];

        xr_k      = xk(1);
        yr_k      = xk(2);
        thetar_k  = xk(3);
        theta_t_k = xk(4);

        % robot COM position
        p_r_k = [xr_k; yr_k];

        % hitch position from robot COM
        p_h_k = [xr_k - d*cos(thetar_k);
                 yr_k - d*sin(thetar_k)];

        % trailer COM position from hitch
        p_tr_k = p_h_k - [Lt*cos(theta_t_k);
                          Lt*sin(theta_t_k)];

        % robot tracking
        if isfield(ref_to,'p_r')
            p_r_ref_k = ref_to.p_r(:,k);
            e_r_k     = p_r_k - p_r_ref_k;
            obj       = obj + e_r_k.' * Qr * e_r_k;
        end

        % trailer tracking
        if isfield(ref_to,'p_tr')
            p_tr_ref_k = ref_to.p_tr(:,k);
            e_tr_k     = p_tr_k - p_tr_ref_k;
            obj        = obj + e_tr_k.' * Qt * e_tr_k;
        end

        % articulation penalty
        phi_k = theta_t_k - thetar_k;
        obj   = obj + w_phi * (phi_k^2);

        % control magnitude
        obj = obj + uk.' * R * uk;

        % control smoothness
        if k > 1
            ukm1 = U(:,k-1);
            duk  = uk - ukm1;
            obj  = obj + duk.' * S * duk;
        end
    end

    %% terminal cost
    xN = X(:,N+1);
    xr_N      = xN(1);
    yr_N      = xN(2);
    thetar_N  = xN(3);
    theta_t_N = xN(4);

    p_r_N = [xr_N; yr_N];
    p_h_N = [xr_N - d*cos(thetar_N);
             yr_N - d*sin(thetar_N)];
    p_tr_N = p_h_N - [Lt*cos(theta_t_N);
                      Lt*sin(theta_t_N)];

    if isfield(ref_to,'p_r')
        p_r_ref_N = ref_to.p_r(:,N+1);
        e_r_N     = p_r_N - p_r_ref_N;
        obj       = obj + e_r_N.' * Qr_f * e_r_N;
    end

    if isfield(ref_to,'p_tr')
        p_tr_ref_N = ref_to.p_tr(:,N+1);
        e_tr_N     = p_tr_N - p_tr_ref_N;
        obj        = obj + e_tr_N.' * Qt_f * e_tr_N;
    end

    % terminal articulation penalty
    phi_N = theta_t_N - thetar_N;
    obj   = obj + w_phi_f * (phi_N^2);

    % terminal yaw alignment
    psi_ref_N = yaw_ref(N+1);
    e_yaw_N   = atan2( sin(thetar_N - psi_ref_N), ...
                       cos(thetar_N - psi_ref_N) );
    obj       = obj + w_yaw_f * (e_yaw_N^2);

    %% pack decision variables
    OPT_vars = [reshape(X, nx*(N+1), 1);
                reshape(U, nu*N,      1)];

    %% build NLP
    nlp = struct;
    nlp.x = OPT_vars;
    nlp.f = obj;
    nlp.g = g;   % only dynamics + initial condition

    opts = struct;
    opts.ipopt.print_level = 3;
    opts.print_time        = 0;
    opts.ipopt.acceptable_tol            = 1e-6;
    opts.ipopt.acceptable_obj_change_tol = 1e-4;

    solver = nlpsol('solver','ipopt',nlp,opts);

    %% equality constraints: g == 0
    ng  = numel(g);
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

    %% variable bounds
    nX = nx*(N+1);
    nU = nu*N;

    lbx = -inf(size(OPT_vars));
    ubx =  inf(size(OPT_vars));

    lbx(1:nX) = repmat(x_min, N+1, 1);
    ubx(1:nX) = repmat(x_max, N+1, 1);

    lbx(nX+1:nX+nU) = repmat(u_min, N, 1);
    ubx(nX+1:nX+nU) = repmat(u_max, N, 1);

    %% initial guess
    X_init = repmat(x0_kin(:), 1, N+1);
    U_init = zeros(nu, N);

    x0_guess = [reshape(X_init, nX, 1);
                reshape(U_init, nU, 1)];

    %% solve NLP
    sol = solver('x0',  x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    % unpack solution
    X_opt = reshape(sol_vec(1:nX), nx, N+1);
    U_opt = reshape(sol_vec(nX+1:end), nu, N);

    x_sol = X_opt.';   % (N+1) x 4
    u_sol = U_opt.';   % N x 2

    % debug
    to_dbg = struct;
    to_dbg.obj     = full(sol.f);
    to_dbg.sol_vec = sol_vec;
    to_dbg.stats   = solver.stats();
end
