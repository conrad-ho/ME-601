function [F_drive, tau_r, qp_dbg] = task_space_qp_controller_proj(x, t, params, yd, ydot_d, yddot_ff)
% task-space QP with projected dynamics
% do not explicitly optimize a, lambda; only optimize u = [F_drive; tau_r]
%
% dependency:
%   dyn = towing_dynamic(x, params, 'task')

    dyn = towing_dynamic(x, params, 'task');
  
    y    = dyn.y;        % 2x1 hitch position
    ydot = dyn.ydot;     % 2x1 hitch velocity
    A_y  = dyn.A_y;      % 2x2, yddot = A_y*u + b_y
    b_y  = dyn.b_y;      % 2x1

    % default trajectory: hold current y
     % 默认轨迹：hold
    if nargin < 4 || isempty(yd),      yd      = y;          end
    if nargin < 5 || isempty(ydot_d),  ydot_d  = [0;0];      end
    if nargin < 6 || isempty(yddot_ff),yddot_ff= [0;0];      end

    % task-space PD 生成期望加速度
    Kp = diag([10, 10]);
    Kd = diag([5, 5]);
    
    yddot_des = yddot_ff ...
                - Kd * (ydot - ydot_d) ...
                - Kp * (y - yd);

    % =========================================================
    % task-space affine model:
    %   yddot(u) = A_y*u + b_y
    % error:
    %   e(u) = (A_y*u + b_y) - yddot_des
    % =========================================================

    Wy = eye(2);
    Wu = diag([1e-3, 1e-1]);   % same as your original

    % QP objective: 0.5*||e||_Wy^2 + 0.5*u'*Wu*u
    Hqp = A_y.' * Wy * A_y + Wu;               % 2x2
    fqp = A_y.' * Wy * (b_y - yddot_des);      % 2x1
    % =========================================================
    % bounds on u
    % =========================================================
    if isfield(params,'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params,'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end

    lb = [-Fmax; -Taumax];
    ub = [ Fmax;  Taumax];
    % =========================================================
     % solve 2D QP with CasADi / qpOASES
    % =========================================================
    import casadi.*

    u_sym = SX.sym('u', 2, 1);
    J_expr = 0.5 * u_sym.' * Hqp * u_sym + fqp.' * u_sym;
    
    qp = struct;
    qp.x = u_sym;
    qp.f = J_expr;
    qp.g = [];  % 无等式/不等式

    opt = struct;
    opt.printLevel = 'low';
    opt.print_time = false;

    evalc("solver = qpsol('qp_solver', 'qpoases', qp, opt);");

    try
        [~, sol] = evalc("solver('x0', zeros(2,1), 'lbx', lb, 'ubx', ub, 'lbg', [], 'ubg', []);");
        u_opt = full(sol.x);
        exitflag = 1;
    catch err
        warning('casadi_qp_failed:solver', 'casadi qp failed: %s', err.message);        u_opt = zeros(2,1);
        exitflag = -1;
    end

    F_drive = u_opt(1);
    tau_r   = u_opt(2);

    if nargout > 2
        a_opt = dyn.a(u_opt);   

        qp_dbg = struct;
        qp_dbg.u_opt   = u_opt;
        qp_dbg.exitflag= exitflag;
        qp_dbg.H       = Hqp;
        qp_dbg.f       = fqp;
        qp_dbg.A_y     = A_y;
        qp_dbg.b_y     = b_y;

        % projected dynamics terms (available in 'task' mode)
        qp_dbg.P_u       = dyn.P_u;
        qp_dbg.p0        = dyn.p0;
        qp_dbg.K         = dyn.K;
        qp_dbg.a_opt     = a_opt;
        
        % residual monitors
        qp_dbg.res_acc   = dyn.res_acc(u_opt);
        qp_dbg.res_vel   = dyn.res_vel();
        qp_dbg.condK     = dyn.condK();
        %{
        % keep your old fields for compatibility
        lambda_opt = zeros(3,1);
        qp_dbg.lambda_opt = lambda_opt;
        qp_dbg.z_opt = [a_opt; u_opt; lambda_opt];  % 11x1 dummy
        qp_dbg.Aeq   = [];
        qp_dbg.beq   = [];
        qp_dbg.A     = [];
        qp_dbg.b     = [];
        %}

    end
end
