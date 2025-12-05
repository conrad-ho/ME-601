function [F_drive, tau_r, qp_dbg] = task_space_qp_controller_proj(x, t, params, yd, ydot_d, yddot_ff)
% task-space QP with projected dynamics
% 再也不显式优化 a, lambda，只优化 u = [F_drive; tau_r]

    d = params.d;

    % === dynamics matrices ===
    dyn   = towing_dynamics_mats(x, params);
    M     = dyn.M;        % 6x6
    B     = dyn.B;        % 6x2
    J     = dyn.J;        % 3x6 (constraints)
    dotJ  = dyn.dotJ;     % 3x6
    v     = dyn.v;        % 6x1
    J_y   = dyn.J_y;      % 2x6 (hitch task)
    Jdot_y_v = dyn.Jdot_y_v;  % 2x1 = Jdot_y * v
    Hgen  = dyn.H;        % 6x1 = H(q,qdot)

    % 当前 hitch 输出
    xr = x(1); yr = x(2); thetar = x(3);
    y = [xr - d*cos(thetar);
         yr - d*sin(thetar)];
    ydot = J_y * v;

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
    % 一步：构造投影动力学  a = P*u + p0
    % 来自：
    %   M a + J^T lambda = B u - H
    %   J a + dotJ v = 0
    % 解出 a(u)
    % =========================================================

    invM = M \ eye(size(M));      % 不要直接 inv(M)

    K = J * invM * J.';           % 3x3, SPD
    % 注意依然用线性求解，不要真 inv(K)
    % 下面会用  K \ something

    % 拆成 Q = B u - H，然后代入你截图里的公式
    % a = invM*(B*u - H) - invM*J' * K^{-1}(J*invM*(B*u - H) + dotJ*v)
    % 线性化：a = P*u + p0
    JinvM  = J * invM;            % 3x6
    invMJt = invM * J.';          % 6x3

    % 与 u 有关的部分
    P_u = invM * B ...
        - invMJt * (K \ (JinvM * B));   % 6x2

    % 常数项 p0
    term0 = - JinvM * Hgen + dotJ * v;  % 3x1 里面的 ( -J*invM*H + dotJ*v )
    p0 = - invM * Hgen ...
         - invMJt * (K \ term0);        % 6x1

    % 现在 a = P_u * u + p0

    % =========================================================
    % 二步：task-space 动力学  yddot = A_y*u + b_y
    % =========================================================
    A_y = J_y * P_u;                 % 2x2
    b_y = J_y * p0 + Jdot_y_v;       % 2x1

    % 误差 e = (A_y*u + b_y) - yddot_des
    Wy = eye(2);
    Wu = 1e-3 * eye(2);              % 对 u 做一点正则

    % 成本 0.5 * ||Wy^(1/2) e||^2 + 0.5 u'Wu u
    H = A_y.' * Wy * A_y + Wu;       % 2x2
    f = A_y.' * Wy * (b_y - yddot_des);  % 2x1

    % =========================================================
    % 约束：只有 u 的上下界
    % =========================================================
    if isfield(params,'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params,'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end

    lb = [-Fmax; -Taumax];
    ub = [ Fmax;  Taumax];

    % =========================================================
    % 用 CasADi / qpOASES 解一个 2 维 QP
    % =========================================================
    import casadi.*

    u_sym = SX.sym('u', 2, 1);
    J_expr = 0.5 * u_sym.' * H * u_sym + f.' * u_sym;

    qp.x = u_sym;
    qp.f = J_expr;
    qp.g = [];  % 无等式/不等式
    opt = struct;
    opt.printLevel = 'none';
    opt.print_time = false;
    solver = qpsol('qp_solver', 'qpoases', qp, opt);

    try
        sol = solver('x0', zeros(2,1), ...
                     'lbx', lb, ...
                     'ubx', ub, ...
                     'lbg', [], 'ubg', []);
        u_opt = full(sol.x);
        exitflag = 1;
    catch err
        warning('casadi_qp_failed:solver', 'casadi qp failed: %s', err.message);        u_opt = zeros(2,1);
        exitflag = -1;
    end

    F_drive = u_opt(1);
    tau_r   = u_opt(2);

    if nargout > 2
         a_opt = P_u * u_opt + p0;  
         lambda_opt = zeros(3,1);
        qp_dbg.u_opt   = u_opt;
        qp_dbg.exitflag= exitflag;
        qp_dbg.H       = H;
        qp_dbg.f       = f;
        qp_dbg.A_y     = A_y;
        qp_dbg.b_y     = b_y;
        qp_dbg.P_u     = P_u;
        qp_dbg.p0      = p0;
        %为了
        qp_dbg.a_opt=a_opt ;
        qp_dbg.u_opt=u_opt;
        qp_dbg.lambda_opt=lambda_opt;
        qp_dbg.exitflag=exitflag;

        qp_dbg.z_opt = [a_opt; u_opt; lambda_opt];  % 11x1 dummy
        qp_dbg.Aeq   = [];  % proj 版没有显式等式约束
        qp_dbg.beq   = [];
        qp_dbg.A     = [];
        qp_dbg.b     = [];
    end
end
