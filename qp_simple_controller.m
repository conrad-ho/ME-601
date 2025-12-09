function [F_drive, tau_r, qp_dbg] = qp_simple_controller(x, t, params)
%QP_SIMPLE_CONTROLLER with sample-and-hold
%   - Builds u_nom from PD on trailer COM
%   - Projects onto box [u_min, u_max]
%   - Only recomputes u every dt_qp seconds (sample-and-hold)
%
%   x      : 12x1 state
%   t      : current time
%   params :
%      .ref_to.t, .ref_to.p_r, .ref_to.p_tr
%      optional:
%        .Kp_tr, .Kd_tr
%        .dt_qp        (controller update period, default 0.02)
%        .qp_verbose   (logical, print debug info)

    % ---------- persistent for sample-and-hold & debug ----------
    persistent call_count last_t last_u

    if isempty(call_count)
        call_count = 0;
    end
    if isempty(last_t)
        last_t = -inf;
    end
    if isempty(last_u)
        last_u = [0; 0];
    end

    call_count = call_count + 1;

    qp_verbose = isfield(params,'qp_verbose') && params.qp_verbose;
    if qp_verbose && mod(call_count, 5000) == 0
        fprintf('[QP] call #%d at t = %.6f (before recompute test)\n', ...
                call_count, t);
    end
    % controller update period
    if isfield(params, 'dt_qp')
        dt_qp = params.dt_qp;
    else
        dt_qp = 0.02;   % 50 Hz 默认采样
    end

    % ---------- basic sanity ----------
    if numel(x) ~= 12
        error('qp_simple_controller:bad_x_dim', ...
              'state x must be 12x1, got size %s', mat2str(size(x)));
    end

    if ~isfield(params, 'ref_to')
        error('qp_simple_controller:missing_ref', ...
              'params.ref_to is required (fields: t, p_r, p_tr).');
    end

    ref_to = params.ref_to;
    if ~isfield(ref_to,'t') || ~isfield(ref_to,'p_r') || ~isfield(ref_to,'p_tr')
        error('qp_simple_controller:bad_ref_struct', ...
              'ref_to must contain t, p_r, p_tr.');
    end

    % ---------- 强制输入范围，小一点 ----------
    u_min = [-50; -5];
    u_max = [ 50;  5];

    % 默认不每次都重算 QP，先假设用上一次的 u
    recompute = false;
    if t >= last_t + 0.5*dt_qp   % 留一点余量防止浮点误差
        recompute = true;
    end

    if ~recompute
        % 直接用上次的 u，减少内层循环的负担
        u = last_u;
        F_drive = u(1);
        tau_r   = u(2);

        % debug 里仍然把误差信息留空或简单标记一下
        qp_dbg = struct();
        qp_dbg.u_nom   = [];   % not recomputed
        qp_dbg.u_sat   = u;
        qp_dbg.e_r     = [];
        qp_dbg.e_tr    = [];
        qp_dbg.p_r     = [x(1); x(2)];
        qp_dbg.p_r_d   = [];
        qp_dbg.p_tr    = [x(7); x(8)];
        qp_dbg.p_tr_d  = [];
        qp_dbg.t_query = [];
        return;
    end

    % ========= 只有在需要重算时才走下面的 heavy 部分 =========

    if qp_verbose && (call_count <= 5 || mod(call_count,5000)==0)
        fprintf('[QP] recompute at call #%d, t = %.6f\n', call_count, t);
    end
    
    % ---------- evaluate reference at time t ----------
    t_grid = ref_to.t(:).';      % 1 x (N+1)
    pr     = ref_to.p_r;         % 2 x (N+1)
    ptr    = ref_to.p_tr;        % 2 x (N+1)

    if size(pr,1) ~= 2 || size(ptr,1) ~= 2
        error('qp_simple_controller:bad_ref_dim', ...
              'ref_to.p_r and p_tr must be 2x(N+1).');
    end

    t0      = t_grid(1);
    t_end   = t_grid(end);
    t_query = min(max(t, t0), t_end);

    prx_d  = interp1(t_grid, pr(1,:), t_query, 'linear');
    pry_d  = interp1(t_grid, pr(2,:), t_query, 'linear');
    trx_d  = interp1(t_grid, ptr(1,:), t_query, 'linear');
    try_d  = interp1(t_grid, ptr(2,:), t_query, 'linear');

    p_r_d  = [prx_d; pry_d];    % 2x1
    p_tr_d = [trx_d; try_d];    % 2x1

    % ---------- current positions ----------
    xr = x(1);  yr = x(2);
    xt = x(7);  yt = x(8);

    p_r  = [xr; yr];
    p_tr = [xt; yt];

    e_r  = p_r  - p_r_d;
    e_tr = p_tr - p_tr_d;

    % ---------- 限制误差幅度，避免 e_tr 太大 ----------
    e_max = 10.0;   % 最多按 10m 的误差来算
    norm_e = norm(e_tr);
    if norm_e > e_max
        e_tr = e_tr * (e_max / max(norm_e, 1e-6));
    end

    % ---------- gains ----------
    if isfield(params, 'Kp_tr')
        Kp_tr = params.Kp_tr;
    else
        Kp_tr = diag([10, 10]);
    end
    if isfield(params, 'Kd_tr')
        Kd_tr = params.Kd_tr;
    else
        Kd_tr = diag([3, 3]);
    end

    % ---------- trailer velocity & error ----------
    vxt = x(10);
    vyt = x(11);
    v_tr   = [vxt; vyt];
    v_tr_d = [0; 0];
    ev_tr  = v_tr - v_tr_d;

    % ---------- nominal control ----------
    f_xy = - Kp_tr * e_tr - Kd_tr * ev_tr;

    theta_r = x(3);
    R = [cos(theta_r);  sin(theta_r)];
    T = [-sin(theta_r); cos(theta_r)];

    F_nom   = dot(f_xy, R);
    tau_nom = dot(f_xy, T);

    u_nom = [F_nom; tau_nom];

    if qp_verbose && (call_count <= 5 || mod(call_count,5000)==0)
        fprintf('      u_nom = [%10.3f, %10.3f]\n', F_nom, tau_nom);
    end

    % ---------- projection ----------
    u = min(max(u_nom, u_min), u_max);

    if qp_verbose && (call_count <= 5 || mod(call_count,5000)==0)
        fprintf('      u_sat = [%10.3f, %10.3f]\n', u(1), u(2));
    end

    if any(~isfinite(u))
        error('qp_simple_controller:NaN_in_u', ...
              'NaN/Inf in u at t=%.6f, u = [%g %g]', t, u(1), u(2));
    end

    % ---------- 更新 sample-and-hold 记忆 ----------
    last_u = u;
    last_t = t;

    % ---------- outputs ----------
    F_drive = u(1);
    tau_r   = u(2);

    % ---------- debug ----------
    qp_dbg = struct();
    qp_dbg.u_nom   = u_nom;
    qp_dbg.u_sat   = u;
    qp_dbg.e_r     = e_r;
    qp_dbg.e_tr    = e_tr;
    qp_dbg.p_r     = p_r;
    qp_dbg.p_r_d   = p_r_d;
    qp_dbg.p_tr    = p_tr;
    qp_dbg.p_tr_d  = p_tr_d;
    qp_dbg.t_query = t_query;
end
