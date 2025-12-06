function [F_drive, tau_r, qp_dbg] = task_space_qp_controller(x, t, params, yd, ydot_d, yddot_ff)
% task-space qp controller skeleton
% maps state x and desired task-space trajectory to [F_drive, tau_r]
% uses quadratic programming: min 0.5 z'H z + f'z subject to linear constraints

    % unpack basic stuff
    d = params.d;

    % get dynamics matrices
    dyn = towing_dynamics_mats(x, params);
    M = dyn.M;
    B = dyn.B;
    J = dyn.J;
    dotJ = dyn.dotJ;
    v = dyn.v;
    J_y = dyn.J_y;
    Jdot_y_v = dyn.Jdot_y_v;

    % current hitch output
    xr = x(1); yr = x(2); thetar = x(3);
    y = [xr - d*cos(thetar);
         yr - d*sin(thetar)];
    ydot = J_y * v;

    % if no desired trajectory is given, hold current hitch position
    if nargin < 4 || isempty(yd)
        yd = y;
    end
    if nargin < 5 || isempty(ydot_d)
        ydot_d = [0; 0];
    end
    if nargin < 6 || isempty(yddot_ff)
        yddot_ff = [0; 0];
    end

    % simple pd in task space to get desired yddot
    Kp = diag([10,10]);   % todo: tune
    Kd = diag([5,5]);     % todo: tune

    yddot_des = yddot_ff ...
                - Kd * (ydot - ydot_d) ...
                - Kp * (y - yd);  % 2x1

    % =========================
    % build qp: min 0.5 z'H z + f'z
    % z = [a(6); u(2); lambda(3)] \in R^11
    % =========================

    nv = 6;      % size of a
    nu = 2;      % size of u
    nl = 3;      % size of lambda
    nz = nv + nu + nl;

    % cost on task-space accel tracking: yddot = J_y * a + Jdot_y_v
    Wy = eye(2);   % todo: tune tracking weight

    % yddot_error = (J_y * a + Jdot_y_v) - yddot_des
    % write yddot_error = C*z - d_vec
    C = [J_y, zeros(2, nu + nl)];      % 2 x 11
    d_vec = yddot_des - Jdot_y_v;      % 2 x 1

    % cost: 0.5 * ||C*z - d_vec||_Wy^2
    H = C' * Wy * C;
    f = -C' * Wy * d_vec;

    % small regularization on u and lambda to keep them small/smooth
    Wu = 1e-3 * eye(nu);
    Wl = 1e-4 * eye(nl);

    % build selection matrices for u and lambda
    Su = [zeros(nu, nv), eye(nu), zeros(nu, nl)];   % pick u from z
    Sl = [zeros(nl, nv + nu), eye(nl)];            % pick lambda from z

    H = H + Su' * Wu * Su + Sl' * Wl * Sl;
    % f can stay the same since center for u, lambda is zero

    % =========================
    % equality constraints Aeq * z = beq
    % 1) dynamics: M*a - B*u - J'*lambda = 0
    % 2) constraint accel: J*a + dotJ*v = 0
    % =========================

    % 1) dynamics constraint
    Aeq_dyn = [M, -B, -J'];             % 6 x 11
    beq_dyn = zeros(6, 1);

    % 2) constraint acceleration
    Aeq_con = [J, zeros(size(J,1), nu + nl)];  % 3 x 11
    beq_con = -dotJ * v;                       % 3 x 1

    Aeq = [Aeq_dyn;
           Aeq_con];             % 9 x 11
    beq = [beq_dyn;
           beq_con];

    % =========================
    % inequality / bounds
    % here we only put bounds on u; others unbounded
    % =========================

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

    lb = -inf(nz, 1);
    ub =  inf(nz, 1);

    % bounds for u in z = [a; u; lambda]
    idx_u = nv + (1:nu);
    
    lb(idx_u(1)) = -Fmax;
    ub(idx_u(1)) =  Fmax;
    lb(idx_u(2)) = -Taumax;
    ub(idx_u(2)) =  Taumax;
    
    % no inequality constraints for now
    A = [];
    b = [];

        % =========================
    % solve qp with casadi (high-level interface)
    % =========================

    import casadi.*

    nz = length(f);  % dimension of z = [a; u; lambda]

    % box bounds
    if isempty(lb)
        lb = -inf(nz,1);
    end
    if isempty(ub)
        ub =  inf(nz,1);
    end

    % ---- CasiDi's decision variable: z_sym ----
    z_sym = SX.sym('z', nz, 1);

    % cost function: J(z) = 0.5 z'H z + f'z
    J_expr = 0.5 * z_sym.' * H * z_sym + f.' * z_sym;

    % constratin Aeq * z = beq
    % in form g(z) = Aeq*z - beq = 0
    if ~isempty(Aeq)
        g_expr = Aeq * z_sym - beq;              % ne x 1
        lbg    = zeros(size(beq));               % =0
        ubg    = zeros(size(beq));               % =0
    else
        g_expr = SX.zeros(0,1);                  % no contraint
        lbg    = [];
        ubg    = [];
    end

    % building QP form in CasiDi's struct
    qp = struct();
    qp.x = z_sym;
    qp.f = J_expr;
    qp.g = g_expr;

    % create QP solver
    opt=struct;
    opts.printLevel='none';
    opts.print_time = false;  
    solver = qpsol('qp_solver', 'qpoases', qp,opts); %%qpoases can be exchange for better performance

    try
        sol = solver( ...
            'x0',  zeros(nz,1), ... % initial guess, could be all zero
            'lbx', lb, ...
            'ubx', ub, ...
            'lbg', lbg, ...
            'ubg', ubg);

        z_opt   = full(sol.x);
        exitflag = 1;
    catch err
        warning('casadi_qp_failed:solver', 'casadi qp failed: %s', err.message);
        z_opt    = [];
        exitflag = -1;
    end

    if exitflag <= 0 || isempty(z_opt)
        % qp failed, you can fall back to something simple (e.g. zero input or fl)
        % here we just set u = 0
        a_opt      = zeros(nv, 1);
        u_opt      = zeros(nu, 1);
        lambda_opt = zeros(nl, 1);
    else
        a_opt      = z_opt(1:nv);
        u_opt      = z_opt(nv+1:nv+nu);
        lambda_opt = z_opt(nv+nu+1:end);
    end

    % unpack control inputs
    F_drive = u_opt(1);
    tau_r   = u_opt(2);

    if nargout > 2
        qp_dbg.z_opt      = z_opt;
        qp_dbg.exitflag   = exitflag;
        qp_dbg.a_opt      = a_opt;
        qp_dbg.u_opt      = u_opt;
        qp_dbg.lambda_opt = lambda_opt;
        qp_dbg.H          = H;
        qp_dbg.f          = f;
        qp_dbg.Aeq        = Aeq;
        qp_dbg.beq        = beq;
        qp_dbg.A          = A;
        qp_dbg.b          = b;
    end
end
