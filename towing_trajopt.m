function [x_sol, u_sol, to_dbg] = towing_trajopt(dt, N, ref_to, x0, params)
% towing_trajopt
% basic multiple shooting trajectory optimization for a robot+trailer system.
% trailer/hitch is the main tracking target, robot has weaker soft constraints.

    import casadi.*

    % parse reference path (hitch path)
    % preferred: ref_to.p_hitch, legacy: ref_to.p_tr
    if isfield(ref_to, 'p_hitch')
        p_hitch_ref = ref_to.p_hitch;   % 2 x (N+1)
    elseif isfield(ref_to, 'p_tr')
        p_hitch_ref = ref_to.p_tr;      % 2 x (N+1)
    else
        error('ref_to must contain field p_hitch or p_tr (2 x (N+1)).');
    end

    %% dimensions and bounds
    nx = 12;    % [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
    nu = 2;     % u = [F; tau]

    % control bounds
    if isfield(params, 'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params, 'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end
    u_min = [-Fmax; -Taumax];
    u_max = [ Fmax;  Taumax];

    % state bounds (keep loose first)
    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% casadi symbols
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    % ================================
    % new unified dynamics interface
    % ================================
    % expects: dyn_full.xdot(u) returns 12x1 state derivative [v; a(u)]
    dyn_full = towing_dynamic(x_sym, params, 'full');
    xdot_sym = dyn_full.xdot(u_sym);

    % wrap into a casadi function for use in multiple shooting
    f_dyn = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});

    %% multiple shooting decision variables
    X = SX.sym('X', nx, N+1);
    U = SX.sym('U', nu, N);

    %% weights
    Qp_tr   = diag([10, 10]);     % hitch tracking (running)
    Qp_tr_f = diag([50, 50]);     % hitch tracking (terminal)

    Qp_r = diag([1, 1]);          % robot soft tracking (unused for now)
    w_r  = 0.1;

    w_phi   = 1.0;                % running hitch angle penalty
    w_phi_f = 5.0;                % terminal hitch angle penalty

    R = diag([0.1, 0.1]);         % control effort
    S = diag([1.0, 1.0]);         % control smoothness (unused for now)

    %% objective and constraints
    obj = 0;
    g   = [];

    % initial condition
    g = [g; X(:,1) - x0(:)];

    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        xkp = X(:,k+1);

        % euler integration constraint
        fk = f_dyn(xk, uk);
        g  = [g; xkp - (xk + dt * fk)];

        % hitch position tracking
        p_ref_k = p_hitch_ref(:,k);
        p_k     = hitch_from_state(xk, params);
        e       = p_k - p_ref_k;
        obj     = obj + e.' * Qp_tr * e;

        % robot soft reference (kept as placeholder, still commented)
        %{
        if isfield(ref_to, 'p_r')
            p_r_ref_k = ref_to.p_r(:,k);
        else
            p_r_ref_k = p_ref_k;
        end
        p_r_k = robot_pos_from_state(xk);
        e_r   = p_r_k - p_r_ref_k;
        obj   = obj + w_r * (e_r.' * Qp_r * e_r);
        %}

        % hitch angle penalty
        phi_k = hitch_angle_from_state(xk);
        obj   = obj + w_phi * (phi_k^2);

        % control effort
        obj = obj + uk.' * R * uk;

        % control smoothness (optional)
        %{
        if k > 1
            duk = uk - U(:,k-1);
            obj = obj + duk.' * S * duk;
        end
        %}
    end

    % terminal cost
    xN = X(:,N+1);

    p_ref_N = p_hitch_ref(:,N+1);
    p_N     = hitch_from_state(xN, params);
    eN      = p_N - p_ref_N;
    obj     = obj + eN.' * Qp_tr_f * eN;

    phi_N   = hitch_angle_from_state(xN);
    obj     = obj + w_phi_f * (phi_N^2);

    %% pack nlp
    OPT_vars = [reshape(X, nx*(N+1), 1);
                reshape(U, nu*N, 1)];

    nlp = struct('x', OPT_vars, 'f', obj, 'g', g);

    opts = struct;
    opts.ipopt.print_level = 3;
    opts.print_time = 0;
    opts.ipopt.acceptable_tol = 1e-6;
    opts.ipopt.acceptable_obj_change_tol = 1e-4;

    solver = nlpsol('solver', 'ipopt', nlp, opts);

    %% bounds on constraints (all equalities)
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
    X_init  = repmat(x0(:), 1, N+1);
    U_init  = zeros(nu, N);
    x0_guess = [reshape(X_init, nX, 1);
                reshape(U_init, nU, 1)];

    %% solve
    sol = solver('x0', x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    X_opt = reshape(sol_vec(1:nX), nx, N+1);
    U_opt = reshape(sol_vec(nX+1:end), nu, N);

    x_sol = X_opt.';   % (N+1) x nx
    u_sol = U_opt.';   % N x nu

    to_dbg = struct;
    to_dbg.obj     = full(sol.f);
    to_dbg.sol_vec = sol_vec;
    to_dbg.stats   = solver.stats();
end

%% ===== helper functions =====

function p_h = hitch_from_state(x, params)
    d  = params.d;
    xr = x(1);  yr = x(2);  thetar = x(3);
    p_h = [xr - d*cos(thetar);
           yr - d*sin(thetar)];
end

function p_r = robot_pos_from_state(x)
    xr = x(1);
    yr = x(2);
    p_r = [xr; yr];
end

function phi = hitch_angle_from_state(x)
    theta_r = x(3);
    theta_t = x(9);
    phi = theta_t - theta_r;
end
