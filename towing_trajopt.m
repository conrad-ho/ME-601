function [x_sol, u_sol, to_dbg] = towing_trajopt(dt, N, ref_to, x0, params)
% towing_trajopt
% multiple shooting TO for robot+trailer with RK4 dynamics
% uses soft hitch-closure penalty (not equality)

    import casadi.*

    %% =========================
    % 0) weights & switches
    % =========================
    w = struct;

    % tracking weights (robot-side hitch -> ref hitch)
    w.Qp_tr   = diag([10, 10]);      % running hitch tracking
    w.Qp_tr_f = diag([50, 50]);      % terminal hitch tracking

    % hitch closure penalty weight
    if isfield(params,'w_hitch'), w.w_hitch = params.w_hitch; else, w.w_hitch = 1e3; end

    % hitch angle penalty (phi = theta_t - theta_r)
    w.w_phi   = 1.0;
    w.w_phi_f = 5.0;

    % control effort
    w.R = diag([0.1, 0.1]);

    % optional control smoothness
    w.use_smooth = false;
    w.S = diag([1.0, 1.0]);

    % optional robot position soft tracking (off by default)
    w.use_robot_track = false;
    w.Qp_r = diag([1, 1]);
    w.w_r  = 0.1;

    %% =========================
    % 1) parse reference path (hitch ref)
    % =========================
    if isfield(ref_to, 'p_hitch')
        p_hitch_ref = ref_to.p_hitch;   % 2 x (N+1)
    elseif isfield(ref_to, 'p_tr')
        p_hitch_ref = ref_to.p_tr;      % 2 x (N+1)
    else
        error('ref_to must contain field p_hitch or p_tr (2 x (N+1)).');
    end

    %% =========================
    % 2) dimensions, bounds
    % =========================
    nx = 12;
    nu = 2;

    if isfield(params, 'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params, 'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end
    u_min = [-Fmax; -Taumax];
    u_max = [ Fmax;  Taumax];

    x_min = -inf(nx,1);
    x_max =  inf(nx,1);

    %% =========================
    % 3) casadi symbols + dynamics
    % =========================
    x_sym = SX.sym('x', nx, 1);
    u_sym = SX.sym('u', nu, 1);

    dyn_full = towing_dynamic(x_sym, params, 'full');
    xdot_sym = dyn_full.xdot(u_sym);
    f_dyn = Function('f_dyn', {x_sym, u_sym}, {xdot_sym});

    %% =========================
    % 4) decision variables
    % =========================
    X = SX.sym('X', nx, N+1);
    U = SX.sym('U', nu, N);

    %% =========================
    % 5) objective + constraints
    % =========================
    obj = 0;
    g   = [];

    % initial condition (hard)
    g = [g; X(:,1) - x0(:)];

    for k = 1:N
        xk  = X(:,k);
        uk  = U(:,k);
        xkp = X(:,k+1);

        % rk4 multiple shooting constraint (hard)
        k1 = f_dyn(xk, uk);
        k2 = f_dyn(xk + (dt/2)*k1, uk);
        k3 = f_dyn(xk + (dt/2)*k2, uk);
        k4 = f_dyn(xk + dt*k3, uk);
        xk_next = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        g = [g; xkp - xk_next];

        % hitch tracking (robot hitch -> reference)
        p_ref_k = p_hitch_ref(:,k);
        p_hr    = hitch_from_state(xk, params);
        e_tr    = p_hr - p_ref_k;
        obj     = obj + e_tr.' * w.Qp_tr * e_tr;

        % soft hitch closure penalty (robot hitch == trailer hitch)
        p_ht  = trailer_hitch_from_state(xk, params);
        e_h   = p_hr - p_ht;
        obj   = obj + w.w_hitch * (e_h.'*e_h);

        % optional robot position track (off by default)
        if w.use_robot_track
            if isfield(ref_to, 'p_r')
                p_r_ref_k = ref_to.p_r(:,k);
            else
                p_r_ref_k = p_ref_k;
            end
            p_r_k = robot_pos_from_state(xk);
            e_r   = p_r_k - p_r_ref_k;
            obj   = obj + w.w_r * (e_r.' * w.Qp_r * e_r);
        end

        % hitch angle penalty (wrapped)
        phi_raw = xk(9) - xk(3);
        phi_k   = atan2(sin(phi_raw), cos(phi_raw));
        obj     = obj + w.w_phi * (phi_k^2);

        % control effort
        obj = obj + uk.' * w.R * uk;

        % optional smoothness
        if w.use_smooth && k > 1
            duk = uk - U(:,k-1);
            obj = obj + duk.' * w.S * duk;
        end
    end

    % terminal cost
    xN = X(:,N+1);

    p_ref_N = p_hitch_ref(:,N+1);
    p_hrN   = hitch_from_state(xN, params);
    eN      = p_hrN - p_ref_N;
    obj     = obj + eN.' * w.Qp_tr_f * eN;

    phi_rawN = xN(9) - xN(3);
    phi_N    = atan2(sin(phi_rawN), cos(phi_rawN));
    obj      = obj + w.w_phi_f * (phi_N^2);

    p_htN = trailer_hitch_from_state(xN, params);
    e_hN  = p_hrN - p_htN;
    obj   = obj + w.w_hitch * (e_hN.'*e_hN);

    %% =========================
    % 6) pack NLP
    % =========================
    OPT_vars = [reshape(X, nx*(N+1), 1);
                reshape(U, nu*N, 1)];
    nlp = struct('x', OPT_vars, 'f', obj, 'g', g);

    opts = struct;
    opts.ipopt.print_level = 5;
    opts.print_time = 0;
    opts.ipopt.max_iter = 5000;
    opts.ipopt.mu_strategy = 'adaptive';
    opts.ipopt.sb = 'yes';

    solver = nlpsol('solver', 'ipopt', nlp, opts);

    %% =========================
    % 7) bounds
    % =========================
    ng  = numel(g);
    lbg = zeros(ng,1);
    ubg = zeros(ng,1);

    nX = nx*(N+1);
    nU = nu*N;

    lbx = -inf(size(OPT_vars));
    ubx =  inf(size(OPT_vars));

    lbx(1:nX) = repmat(x_min, N+1, 1);
    ubx(1:nX) = repmat(x_max, N+1, 1);

    lbx(nX+1:nX+nU) = repmat(u_min, N, 1);
    ubx(nX+1:nX+nU) = repmat(u_max, N, 1);

    %% =========================
    % 8) initial guess
    % =========================
    X_init = repmat(x0(:), 1, N+1);
    U_init = zeros(nu, N);

    % seed robot position from ref hitch (keeps angles/velocities from x0)
    try
        for k = 1:(N+1)
            th = X_init(3,k);
            X_init(1,k) = p_hitch_ref(1,k) + params.d*cos(th);
            X_init(2,k) = p_hitch_ref(2,k) + params.d*sin(th);
        end
    catch
    end

    x0_guess = [reshape(X_init, nX, 1);
                reshape(U_init, nU, 1)];

    %% =========================
    % 9) solve
    % =========================
    sol = solver('x0', x0_guess, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);

    sol_vec = full(sol.x);

    X_opt = reshape(sol_vec(1:nX), nx, N+1);
    U_opt = reshape(sol_vec(nX+1:end), nu, N);

    x_sol = X_opt.';   % (N+1) x nx
    u_sol = U_opt.';   % N x nu

    %% =========================
    % 10) debug + cost monitors (post-eval, does not change solution)
    % =========================
    to_dbg = struct;
    to_dbg.obj     = full(sol.f);
    to_dbg.sol_vec = sol_vec;
    to_dbg.stats   = solver.stats();
    to_dbg.w       = w;
    to_dbg.dt      = dt;
    to_dbg.p_hitch_ref = p_hitch_ref;

    % per-step monitors
    J_pos = zeros(N,1);
    J_vel = zeros(N,1);
    J_u   = zeros(N,1);
    J_yaw = zeros(N,1);
    J_h   = zeros(N,1);
    J_sm  = zeros(N,1);

    % hitch positions for finite-diff ydot
    p_hr_all = zeros(2, N+1);
    for k = 1:(N+1)
        p_hr_all(:,k) = hitch_from_state(x_sol(k,:).', params);
    end

    % velocity weight for monitor only (you can change later)
    Qv = w.Qp_tr;

    for k = 1:N
        xk = x_sol(k,:).';
        uk = u_sol(k,:).';

        % pos tracking term
        p_ref_k = p_hitch_ref(:,k);
        p_hr    = p_hr_all(:,k);
        e_tr    = p_hr - p_ref_k;
        J_pos(k)= full(e_tr.' * w.Qp_tr * e_tr);

        % vel monitor (finite difference)
        ydot_k     = (p_hr_all(:,k+1) - p_hr_all(:,k)) / dt;
        ydot_ref_k = (p_hitch_ref(:,k+1) - p_hitch_ref(:,k)) / dt;
        e_v        = ydot_k - ydot_ref_k;
        J_vel(k)   = full(e_v.' * Qv * e_v);

        % yaw/phi penalty
        phi_raw = xk(9) - xk(3);
        phi_k   = atan2(sin(phi_raw), cos(phi_raw));
        J_yaw(k)= full(w.w_phi * (phi_k^2));

        % control effort
        J_u(k)  = full(uk.' * w.R * uk);

        % hitch closure penalty
        p_ht = trailer_hitch_from_state(xk, params);
        e_h  = p_hr - p_ht;
        J_h(k)= full(w.w_hitch * (e_h.'*e_h));

        % optional smoothness
        if w.use_smooth && k > 1
            duk   = uk - u_sol(k-1,:).';
            J_sm(k)= full(duk.' * w.S * duk);
        end
    end

    % terminal terms
    xN_num = x_sol(N+1,:).';
    p_hrN  = hitch_from_state(xN_num, params);
    eN     = p_hrN - p_hitch_ref(:,N+1);
    J_pos_f = full(eN.' * w.Qp_tr_f * eN);

    phi_rawN = xN_num(9) - xN_num(3);
    phi_N    = atan2(sin(phi_rawN), cos(phi_rawN));
    J_yaw_f  = full(w.w_phi_f * (phi_N^2));

    p_htN = trailer_hitch_from_state(xN_num, params);
    e_hN  = p_hrN - p_htN;
    J_h_f = full(w.w_hitch * (e_hN.'*e_hN));

    to_dbg.cost = struct;
    to_dbg.cost.J_pos = J_pos;
    to_dbg.cost.J_vel = J_vel; % monitor only unless you add it to obj
    to_dbg.cost.J_u   = J_u;
    to_dbg.cost.J_yaw = J_yaw;
    to_dbg.cost.J_h   = J_h;
    to_dbg.cost.J_sm  = J_sm;
    to_dbg.cost.term  = struct('J_pos_f',J_pos_f,'J_yaw_f',J_yaw_f,'J_h_f',J_h_f);

    to_dbg.cost.sum = struct;
    to_dbg.cost.sum.pos = sum(J_pos) + J_pos_f;
    to_dbg.cost.sum.vel = sum(J_vel);
    to_dbg.cost.sum.u   = sum(J_u);
    to_dbg.cost.sum.yaw = sum(J_yaw) + J_yaw_f;
    to_dbg.cost.sum.h   = sum(J_h) + J_h_f;
    to_dbg.cost.sum.sm  = sum(J_sm);
    to_dbg.cost.sum.total_like = to_dbg.cost.sum.pos + to_dbg.cost.sum.u + to_dbg.cost.sum.yaw + to_dbg.cost.sum.h + to_dbg.cost.sum.sm;

    % closure sanity check (meters)
    try
        xr = x_sol(:,1); yr = x_sol(:,2); thr = x_sol(:,3);
        xt = x_sol(:,7); yt = x_sol(:,8); tht = x_sol(:,9);
        p_hr = [xr - params.d*cos(thr), yr - params.d*sin(thr)];
        p_ht = [xt + (params.Lt/2)*cos(tht), yt + (params.Lt/2)*sin(tht)];
        e_h  = sqrt(sum((p_hr - p_ht).^2,2));
        to_dbg.hitch_err_max = max(e_h);
        to_dbg.hitch_err_rms = rms(e_h);
        fprintf('to hitch closure (soft): max=%.3e m, rms=%.3e m\n', to_dbg.hitch_err_max, to_dbg.hitch_err_rms);
    catch
    end
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

function p_ht = trailer_hitch_from_state(x, params)
    Lt = params.Lt;
    xt = x(7);  yt = x(8);  thetat = x(9);
    p_ht = [xt + (Lt/2)*cos(thetat);
            yt + (Lt/2)*sin(thetat)];
end
