function [x_sol, u_sol, to_dbg] = towing_trajopt(dt, N, ref_to, x0, params)
% towing_trajopt_rigid_penalty
%   multiple shooting TO for robot+trailer.
%   uses towing_dynamics_full (with baumgarte etc.) + rk4 discretization.
%   rigid hitch constraint is enforced softly via a large penalty on
%   position mismatch between robot-side and trailer-side hitch.
%
% state:
%   x = [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt] (12x1)
% input:
%   u = [F_drive; tau_r] (2x1)
%
% geometry (robot on the right, trailer on the left, both facing +x):
%   hitch from robot com:
%       p_hr = [xr - d*cos(thetar);
%               yr - d*sin(thetar)]
%   hitch from trailer com:
%       p_ht = [xt + Lt*cos(thetat);
%               yt + Lt*sin(thetat)]
%   if your trailer hitch is actually at -Lt, change +Lt to -Lt below.
%
% inputs:
%   dt      : time step
%   N       : horizon length (# of control steps, states are N+1)
%   ref_to  : struct (all fields optional)
%               .t    : 1 x (N+1) time grid
%               .p_r  : 2 x (N+1) robot com ref path
%               .p_tr : 2 x (N+1) trailer com ref path
%   x0      : initial state (12x1)
%   params  : struct with fields
%               .d, .Lt, Fmax, Taumax, ...
%
% outputs:
%   x_sol   : (N+1) x 12 optimal state trajectory
%   u_sol   : N x 2    optimal input trajectory
%   to_dbg  : struct with debug info

    import casadi.*

    %% dimensions
    nx = 12;
    nu = 2;

    %% control bounds
    if isfield(params,'Fmax')
        Fmax = params.Fmax;
    else
        Fmax = 1e3;
    end
    if isfield(params,'Taumax')
        Taumax = params.Taumax;
    else
        Taumax = 1e3;
    end
    u_min = [-Fmax; -Taumax];
    u_max = [ Fmax;  Taumax];

    %% state bounds (can be tightened later)
    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% continuous-time dynamics: towing_dynamics_full
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    xdot_sym = towing_dynamics_full(x_sym, u_sym, params);
    f_dyn    = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});

    %% rk4 one-step map: x_{k+1} = F_rk4(x_k,u_k)
    xk_sym = SX.sym('xk', nx, 1);
    uk_sym = SX.sym('uk', nu, 1);

    k1 = f_dyn(xk_sym,          uk_sym);
    k2 = f_dyn(xk_sym + 0.5*dt*k1, uk_sym);
    k3 = f_dyn(xk_sym + 0.5*dt*k2, uk_sym);
    k4 = f_dyn(xk_sym + dt*k3,     uk_sym);

    xk_rk4 = xk_sym + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    F_rk4  = Function('F_rk4', {xk_sym, uk_sym}, {xk_rk4});

    %% hitch position mismatch phi_pos(x) = p_hr - p_ht (position level)
    d  = params.d;
    Lt = params.Lt;

    xr = x_sym(1);  yr = x_sym(2);  thetar  = x_sym(3);
    xt = x_sym(7);  yt = x_sym(8);  thetat  = x_sym(9);

    % hitch from robot com
    p_hr_sym = [xr - d*cos(thetar);
                yr - d*sin(thetar)];

    % hitch from trailer com
    % if your trailer hitch is at -Lt instead of +Lt, change +Lt to -Lt here
    p_ht_sym = [xt + Lt*cos(thetat);
                yt + Lt*sin(thetat)];

    phi_pos_sym = p_hr_sym - p_ht_sym;   % 2x1

    hitch_pos_mismatch = Function('phi_pos', {x_sym}, {phi_pos_sym});

    %% define yaw_ref from robot reference path if available
    yaw_ref = zeros(1, N+1);
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

    %% multiple-shooting variables
    X = SX.sym('X', nx, N+1);   % states
    U = SX.sym('U', nu, N);     % controls

    %% cost weights
    % robot / hitch tracking via robot com
    Qh    = diag([30, 30]);
    Qh_f  = diag([100, 100]);

    % trailer com tracking
    Qtr   = 0.3 * Qh;
    Qtr_f = 0.4 * Qh_f;

    % articulation angle penalty
    w_phi   = 50;
    w_phi_f = 200;

    % hitch position mismatch penalty (rigid constraint soft enforcement)
    w_phi_pos   = 1e4;   % running
    w_phi_pos_f = 1e4;   % terminal (can be larger if needed)

    % control regularization
    R  = diag([0.1, 0.1]);   % magnitude
    S  = diag([0.05, 0.05]); % smoothness

    % yaw regularization
    w_wr    = 0.01;
    w_yaw_f = 0.10;

    %% objective and constraints
    obj = 0;
    g   = [];

    % initial condition
    g = [g; X(:,1) - x0(:)];

    %% main loop
    for k = 1:N
        xk = X(:,k);
        uk = U(:,k);
        xkp = X(:,k+1);

        % ----- discrete dynamics (rk4) -----
        xk_next = F_rk4(xk, uk);
        g = [g; xkp - xk_next];

        % ----- robot com tracking (if ref_to.p_r exists) -----
        if isfield(ref_to,'p_r')
            p_r_ref_k = ref_to.p_r(:,k);
            p_r_k     = [xk(1); xk(2)];
            e_r_k     = p_r_k - p_r_ref_k;
            obj       = obj + e_r_k.' * Qh * e_r_k;
        end

        % ----- trailer com tracking (if ref_to.p_tr exists) -----
        if isfield(ref_to,'p_tr')
            p_tr_ref_k = ref_to.p_tr(:,k);
            p_tr_k     = [xk(7); xk(8)];
            e_tr_k     = p_tr_k - p_tr_ref_k;
            obj        = obj + e_tr_k.' * Qtr * e_tr_k;
        end

        % ----- articulation angle penalty (phi = theta_t - theta_r) -----
        theta_r_k = xk(3);
        theta_t_k = xk(9);
        phi_k     = theta_t_k - theta_r_k;
        obj       = obj + w_phi * (phi_k^2);

        % ----- hitch position mismatch penalty (soft rigid constraint) -----
        phi_pos_k = hitch_pos_mismatch(xk);  % 2x1
        obj       = obj + w_phi_pos * (phi_pos_k.' * phi_pos_k);

        % ----- yaw-rate regularization -----
        wr_k = xk(6);
        obj  = obj + w_wr * (wr_k^2);

        % ----- control magnitude -----
        obj  = obj + uk.' * R * uk;

        % ----- control smoothness -----
        if k > 1
            ukm1 = U(:,k-1);
            duk  = uk - ukm1;
            obj  = obj + duk.' * S * duk;
        end
    end

    %% terminal cost
    xN = X(:,N+1);

    % terminal robot com tracking
    if isfield(ref_to,'p_r')
        p_r_ref_N = ref_to.p_r(:,N+1);
        p_r_N     = [xN(1); xN(2)];
        e_r_N     = p_r_N - p_r_ref_N;
        obj       = obj + e_r_N.' * Qh_f * e_r_N;
    end

    % terminal trailer com tracking
    if isfield(ref_to,'p_tr')
        p_tr_ref_N = ref_to.p_tr(:,N+1);
        p_tr_N     = [xN(7); xN(8)];
        e_tr_N     = p_tr_N - p_tr_ref_N;
        obj        = obj + e_tr_N.' * Qtr_f * e_tr_N;
    end

    % terminal articulation angle penalty
    theta_r_N = xN(3);
    theta_t_N = xN(9);
    phi_N     = theta_t_N - theta_r_N;
    obj       = obj + w_phi_f * (phi_N^2);

    % terminal hitch position mismatch penalty
    phi_pos_N = hitch_pos_mismatch(xN);
    obj       = obj + w_phi_pos_f * (phi_pos_N.' * phi_pos_N);

    % terminal yaw alignment
    theta_r_N = xN(3);
    psi_ref_N = yaw_ref(N+1);
    e_yaw_N   = atan2( sin(theta_r_N - psi_ref_N), ...
                       cos(theta_r_N - psi_ref_N) );
    obj       = obj + w_yaw_f * (e_yaw_N^2);

    %% pack decision variables
    OPT_vars = [reshape(X, nx*(N+1), 1);
                reshape(U, nu*N,      1)];

    %% build nlp
    nlp = struct;
    nlp.x = OPT_vars;
    nlp.f = obj;
    nlp.g = g;    % only dynamics + initial condition (equalities)

    opts = struct;
    opts.ipopt.print_level = 3;
    opts.print_time        = 0;
    opts.ipopt.acceptable_tol            = 1e-6;
    opts.ipopt.acceptable_obj_change_tol = 1e-4;

    solver = nlpsol('solver', 'ipopt', nlp, opts);

    %% constraint bounds: all g are equalities
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
    X_init = repmat(x0(:), 1, N+1);
    U_init = zeros(nu, N);

    x0_guess = [reshape(X_init, nX, 1);
                reshape(U_init, nU, 1)];

    %% solve nlp
    sol = solver('x0',  x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    % unpack solution
    X_opt = reshape(sol_vec(1:nX), nx, N+1);
    U_opt = reshape(sol_vec(nX+1:end), nu, N);

    x_sol = X_opt.';   % (N+1) x 12
    u_sol = U_opt.';   % N x 2

    % debug info
    to_dbg = struct;
    to_dbg.obj     = full(sol.f);
    to_dbg.sol_vec = sol_vec;
    to_dbg.stats   = solver.stats();
end
