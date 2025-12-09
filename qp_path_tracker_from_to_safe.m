function [F_drive, tau_r, qp_dbg] = qp_path_tracker_from_to_safe(x, t, params)
%QP_PATH_TRACKER_FROM_TO_SAFE
%   Very conservative path-tracking controller that only uses TO path:
%     - ref_to.t, ref_to.p_tr, ref_to.p_r
%   and full state x, with strong clipping on errors and inputs.
%
%   x      : 12x1 state
%            [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
%   t      : current time
%   params : struct, must contain:
%              .ref_to.t    [1 x (N+1)]
%              .ref_to.p_tr [2 x (N+1)]
%              .ref_to.p_r  [2 x (N+1)]
%            optional:
%              .u_min (2x1), .u_max (2x1)
%              .Kp_tr (2x2), .Kd_tr (2x2)

    % ---------- sanity ----------
    if numel(x) ~= 12
        error('qp_path_tracker_from_to_safe:bad_x_dim', ...
              'state x must be 12x1, got size %s', mat2str(size(x)));
    end

    if ~isfield(params, 'ref_to')
        error('qp_path_tracker_from_to_safe:missing_ref', ...
              'params.ref_to is required (with fields t, p_tr, p_r).');
    end

    ref_to = params.ref_to;
    if ~isfield(ref_to, 't') || ~isfield(ref_to, 'p_tr') || ~isfield(ref_to, 'p_r')
        error('qp_path_tracker_from_to_safe:bad_ref_struct', ...
              'ref_to must contain t, p_tr, p_r.');
    end

    % ---------- very small input bounds ----------
    if isfield(params, 'u_min')
        u_min = params.u_min(:);
    else
        u_min = [-5; -1];     % 很保守
    end

    if isfield(params, 'u_max')
        u_max = params.u_max(:);
    else
        u_max = [ 5;  1];
    end

    if numel(u_min) ~= 2 || numel(u_max) ~= 2
        error('qp_path_tracker_from_to_safe:bad_bounds', ...
              'u_min and u_max must be 2x1 vectors.');
    end

    % ---------- interpolate reference ----------
    t_grid    = ref_to.t(:).';
    p_tr_ref  = ref_to.p_tr;
    p_r_ref   = ref_to.p_r;

    if size(p_tr_ref,1) ~= 2 || size(p_r_ref,1) ~= 2
        error('qp_path_tracker_from_to_safe:bad_ref_dim', ...
              'ref_to.p_tr and ref_to.p_r must be 2x(N+1).');
    end

    t0    = t_grid(1);
    t_end = t_grid(end);
    t_q   = min(max(t, t0), t_end);

    xtr_d = interp1(t_grid, p_tr_ref(1,:), t_q, 'linear');
    ytr_d = interp1(t_grid, p_tr_ref(2,:), t_q, 'linear');
    xr_d  = interp1(t_grid, p_r_ref(1,:),  t_q, 'linear');
    yr_d  = interp1(t_grid, p_r_ref(2,:),  t_q, 'linear');

    p_tr_d = [xtr_d; ytr_d];
    p_r_d  = [xr_d;  yr_d];

    % ---------- current positions ----------
    xr = x(1);  yr = x(2);
    xt = x(7);  yt = x(8);

    p_r  = [xr; yr];
    p_tr = [xt; yt];

    e_tr = p_tr - p_tr_d;
    e_r  = p_r  - p_r_d; %#ok<NASGU>

    % ---------- clip error magnitude (very small) ----------
    e_max = 1.0;   % 最多按 1m 误差算
    e_norm = norm(e_tr);
    if e_norm > e_max
        e_tr = e_tr * (e_max / max(e_norm, 1e-6));
    end

    % ---------- small gains ----------
    if isfield(params, 'Kp_tr')
        Kp_tr = params.Kp_tr;
    else
        Kp_tr = diag([3, 3]);   % 比之前小很多
    end

    if isfield(params, 'Kd_tr')
        Kd_tr = params.Kd_tr;
    else
        Kd_tr = diag([1, 1]);
    end

    % ---------- trailer velocity ----------
    vxt = x(10);
    vyt = x(11);
    v_tr   = [vxt; vyt];
    v_tr_d = [0; 0];
    ev_tr  = v_tr - v_tr_d;

    % 也把速度误差 clip 一下，防止太大
    v_max = 2.0;
    v_norm = norm(ev_tr);
    if v_norm > v_max
        ev_tr = ev_tr * (v_max / max(v_norm, 1e-6));
    end

    % ---------- virtual planar force ----------
    f_xy = - Kp_tr * e_tr - Kd_tr * ev_tr;

    % 再 clip 一下 f_xy 本身，避免虚拟力太大
    f_max = 20;
    f_norm = norm(f_xy);
    if f_norm > f_max
        f_xy = f_xy * (f_max / max(f_norm, 1e-6));
    end

    % ---------- map to [F_drive; tau_r] ----------
    theta_r = x(3);
    R = [cos(theta_r);  sin(theta_r)];
    T = [-sin(theta_r); cos(theta_r)];

    F_nom   = dot(f_xy, R);
    tau_nom = dot(f_xy, T);
    u_nom   = [F_nom; tau_nom];

    % ---------- "QP": project onto box ----------
    u = min(max(u_nom, u_min), u_max);

    F_drive = u(1);
    tau_r   = u(2);

    % ---------- debug ----------
    qp_dbg = struct();
    qp_dbg.t_query = t_q;
    qp_dbg.p_tr    = p_tr;
    qp_dbg.p_tr_d  = p_tr_d;
    qp_dbg.p_r     = p_r;
    qp_dbg.p_r_d   = p_r_d;
    qp_dbg.e_tr    = e_tr;
    qp_dbg.u_nom   = u_nom;
    qp_dbg.u_sat   = u;
end
