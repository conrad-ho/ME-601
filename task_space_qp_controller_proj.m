function [F_drive, tau_r, qp_dbg] = task_space_qp_controller_proj(x, t, params, yd, ydot_d, yddot_ff)
% task-space QP with projected dynamics (analytic 2D QP, no CasADi)
% we minimize:
%   0.5 * (A_y*u + b_y - yddot_des)' * Wy * (A_y*u + b_y - yddot_des) + 0.5 * u'Wu*u
% subject to:
%   lb <= u <= ub
%
% inputs:
%   x        : 12x1 state
%   t        : time (unused)
%   params   : struct with fields d, Fmax, Taumax, ...
%   yd       : 2x1 desired hitch position
%   ydot_d   : 2x1 desired hitch velocity
%   yddot_ff : 2x1 feedforward hitch acceleration
%
% outputs:
%   F_drive, tau_r : scalar inputs
%   qp_dbg         : debug struct (optional)

    d = params.d;

    % === dynamics matrices ===
    dyn      = towing_dynamics_mats(x, params);
    M        = dyn.M;          % 6x6
    B        = dyn.B;          % 6x2
    J        = dyn.J;          % 3x6 (constraints)
    dotJ     = dyn.dotJ;       % 3x6
    v        = dyn.v;          % 6x1
    J_y      = dyn.J_y;        % 2x6 (hitch task)
    Jdot_y_v = dyn.Jdot_y_v;   % 2x1 = Jdot_y * v
    Hgen     = dyn.H;          % 6x1 = H(q,qdot)

    % current hitch output
    xr = x(1); yr = x(2); thetar = x(3);
    y = [xr - d*cos(thetar);
         yr - d*sin(thetar)];
    ydot = J_y * v;

    % default references (hold)
    if nargin < 4 || isempty(yd),       yd       = y;      end
    if nargin < 5 || isempty(ydot_d),   ydot_d   = [0;0];  end
    if nargin < 6 || isempty(yddot_ff), yddot_ff = [0;0];  end

    % task-space PD desired acceleration
    Kp = diag([10, 10]);
    Kd = diag([5, 5]);

    yddot_des = yddot_ff ...
                - Kd * (ydot - ydot_d) ...
                - Kp * (y   - yd);

    % =========================================================
    % step 1: projected dynamics  a = P_u * u + p0
    % from:
    %   M a + J^T lambda = B u - H
    %   J a + dotJ v = 0
    % =========================================================

    invM = M \ eye(size(M));      % 6x6, numeric solve instead of inv(M)

    K = J * invM * J.';           % 3x3, SPD

    % helper products
    JinvM  = J * invM;            % 3x6
    invMJt = invM * J.';          % 6x3

    % u-dependent part
    P_u = invM * B ...
        - invMJt * (K \ (JinvM * B));        % 6x2

    % constant term p0
    term0 = - JinvM * Hgen + dotJ * v;       % 3x1
    p0 = - invM * Hgen ...
         - invMJt * (K \ term0);             % 6x1

    % so a = P_u * u + p0

    % =========================================================
    % step 2: task-space dynamics  yddot = A_y*u + b_y
    % =========================================================
    A_y = J_y * P_u;                 % 2x2
    b_y = J_y * p0 + Jdot_y_v;       % 2x1

    % error e = (A_y*u + b_y) - yddot_des
    Wy = eye(2);
    Wu = 1e-3 * eye(2);              % small regularization on u

    % quadratic cost in standard form: 0.5 u' H u + f' u
    H = A_y.' * Wy * A_y + Wu;               % 2x2, SPD
    f = A_y.' * Wy * (b_y - yddot_des);      % 2x1

    % =========================================================
    % bounds on u
    % =========================================================
    if isfield(params,'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params,'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end

    lb = [-Fmax; -Taumax];
    ub = [ Fmax;  Taumax];

    % =========================================================
    % step 3: analytic solution of 2D bound QP
    %   min  0.5 u' H u + f' u
    %   s.t. lb <= u <= ub
    %
    % strategy:
    %   - unconstrained optimum u* = -H^{-1} f
    %   - 1D faces (fix u1, optimize u2; fix u2, optimize u1)
    %   - 4 corners
    %   pick the one with minimum cost
    % =========================================================

    % 1) unconstrained optimum
    u_uncon = - H \ f;   % 2x1

    cost = @(u) 0.5*u.'*H*u + f.'*u;

    cand_u   = zeros(2, 0);
    cand_val = [];

    % candidate 1: unconstrained optimum if inside box
    if all(u_uncon >= lb) && all(u_uncon <= ub)
        cand_u   = [cand_u, u_uncon];
        cand_val = [cand_val, cost(u_uncon)];
    end

    % extract entries for convenience
    H11 = H(1,1); H12 = H(1,2);
    H21 = H(2,1); H22 = H(2,2);
    f1  = f(1);   f2  = f(2);

    % 2) 1D faces: fix u1 at bounds, optimize u2
    for u1_fixed = [lb(1), ub(1)]
        if abs(H22) > 1e-12
            u2_star = -(H21*u1_fixed + f2)/H22;
            % clamp to [lb2, ub2]
            u2_star = min(max(u2_star, lb(2)), ub(2));
            u_cand = [u1_fixed; u2_star];
            cand_u   = [cand_u, u_cand];
            cand_val = [cand_val, cost(u_cand)];
        end
    end

    % 3) 1D faces: fix u2 at bounds, optimize u1
    for u2_fixed = [lb(2), ub(2)]
        if abs(H11) > 1e-12
            u1_star = -(H12*u2_fixed + f1)/H11;
            % clamp to [lb1, ub1]
            u1_star = min(max(u1_star, lb(1)), ub(1));
            u_cand = [u1_star; u2_fixed];
            cand_u   = [cand_u, u_cand];
            cand_val = [cand_val, cost(u_cand)];
        end
    end

    % 4) 4 corners
    corners = [lb(1) lb(1) ub(1) ub(1);
               lb(2) ub(2) lb(2) ub(2)];
    for k = 1:4
        u_cand = corners(:,k);
        cand_u   = [cand_u, u_cand];
        cand_val = [cand_val, cost(u_cand)];
    end

    % pick best candidate
    [~, idx_best] = min(cand_val);
    u_opt = cand_u(:, idx_best);

    F_drive = u_opt(1);
    tau_r   = u_opt(2);
    exitflag = 1;

    % =========================================================
    % optional debug info
    % =========================================================
    if nargout > 2
        a_opt       = P_u * u_opt + p0;    % 6x1
        lambda_opt  = zeros(3,1);         % projected version does not solve lambda explicitly

        qp_dbg.u_opt       = u_opt;
        qp_dbg.exitflag    = exitflag;
        qp_dbg.H           = H;
        qp_dbg.f           = f;
        qp_dbg.A_y         = A_y;
        qp_dbg.b_y         = b_y;
        qp_dbg.P_u         = P_u;
        qp_dbg.p0          = p0;
        qp_dbg.a_opt       = a_opt;
        qp_dbg.lambda_opt  = lambda_opt;
        qp_dbg.z_opt       = [a_opt; u_opt; lambda_opt];  % 11x1 dummy
        qp_dbg.Aeq         = [];
        qp_dbg.beq         = [];
        qp_dbg.A           = [];
        qp_dbg.b           = [];
    end
end
