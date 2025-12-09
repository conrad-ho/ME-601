function [F_drive, tau_r, qp_dbg] = qp_path_tracker_from_to(x, t, params)
%QP_PATH_TRACKER_FROM_TO
%   Path-tracking "QP" controller that only uses TO path (ref_to)
%   and full dynamics state x.
%
%   x      : 12x1 state
%            [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
%   t      : current time
%   params : struct, must contain:
%              .ref_to.t    [1 x (N+1)]
%              .ref_to.p_tr [2 x (N+1)]  trailer COM path
%              .ref_to.p_r  [2 x (N+1)]  robot COM path
%            optional:
%              .u_min (2x1), .u_max (2x1)
%              .Kp_tr (2x2), .Kd_tr (2x2)
%
%   Outputs:
%      F_drive : scalar drive force
%      tau_r   : scalar steering torque
%      qp_dbg  : struct (for debug/plotting, can be ignored)

    % ---------- sanity checks ----------
    if numel(x) ~= 12
        error('qp_path_tracker_from_to:bad_x_dim', ...
              'state x must be 12x1, got size %s', mat2str(size(x)));
    end

    if ~isfield(params, 'ref_to')
        error('qp_path_tracker_from_to:missing_ref', ...
              'params.ref_to is required (with fields t, p_tr, p_r).');
    end

    ref_to = params.ref_to;
    if ~isfield(ref_to, 't') || ~isfield(ref_to, 'p_tr') || ~isfield(ref_to, 'p_r')
        error('qp_path_tracker_from_to:bad_ref_struct', ...
              'ref_to must contain t, p_tr, p_r.');
    end

    % ---------- input bounds (safe defaults) ----------
    % 之后你可以根据仿真效果调大/调小
    if isfield(params, 'u_min')
        u_min = params.u_min(:);
    else
        u_min = [-30; -5];     % [F_min; tau_min]
    end

    if isfield(params, 'u_max')
        u_max = params.u_max(:);
    else
        u_max = [ 30;  5];     % [F_max; tau_max]
    end

    if numel(u_min) ~= 2 || numel(u_max) ~= 2
        error('qp_path_tracker_from_to:bad_bounds', ...
              'u_min and u_max must be 2x1 vectors.');
    end

    % ---------- interpolate reference positions at time t ----------
    t_grid = ref_to.t(:).';      % 1 x (N+1)
    p_tr_ref = ref_to.p_tr;      % 2 x (N+1)
    p_r_ref  = ref_to.p_r;       % 2 x (N+1)

    if size(p_tr_ref,1) ~= 2 || size(p_r_ref,1) ~= 2
        error('qp_path_tracker_from_to:bad_ref_dim', ...
              'ref_to.p_tr and ref_to.p_r must be 2x(N+1).');
    end

    % time clamp
    t0    = t_grid(1);
    t_end = t_grid(end);
    t_q   = min(max(t, t0), t_end);

    % trailer desired pos
    xtr_d = interp1(t_grid, p_tr_ref(1,:), t_q, 'linear');
    ytr_d = interp1(t_grid, p_tr_ref(2,:), t_q, 'linear');

    % robot desired pos（现在没直接用到，但放进 dbg 里）
    xr_d  = interp1(t_grid, p_r_ref(1,:),  t_q, 'linear');
    yr_d  = interp1(t_grid, p_r_ref(2,:),  t_q, 'linear');

    p_tr_d = [xtr_d; ytr_d];
    p_r_d  = [xr_d;  yr_d];

    if any(~isfinite(p_tr_d)) || any(~isfinite(p_r_d))
        error('qp_path_tracker_from_to:NaN_in_ref', ...
              'NaN/Inf in interpolated ref at t = %.6f', t_q);
    end

    % ---------- current positions from state ----------
    xr = x(1);  yr = x(2);
    xt = x(7);  yt = x(8);

    p_r  = [xr; yr];
    p_tr = [xt; yt];

    % position errors
    e_tr = p_tr - p_tr_d;
    e_r  = p_r  - p_r_d;   %#ok<NASGU> % 目前没直接用，可以以后加 robot cost

    % ---------- optional: limit error magnitude ----------
    % 避免一开始跑偏太远的时候，误差过大导致 f_xy 巨大 → 数值变 stiff
    e_max = 5.0;          % 最多按 5m 误差来算
    e_norm = norm(e_tr);
    if e_norm > e_max
        e_tr = e_tr * (e_max / max(e_norm, 1e-6));
    end

    % ---------- gains ----------
    if isfield(params, 'Kp_tr')
        Kp_tr = params.Kp_tr;
    else
        Kp_tr = diag([10, 10]);   % position gains
    end

    if isfield(params, 'Kd_tr')
        Kd_tr = params.Kd_tr;
    else
        Kd_tr = diag([3, 3]);     % velocity gains
    end

    if ~isequal(size(Kp_tr), [2 2]) || ~isequal(size(Kd_tr), [2 2])
        error('qp_path_tracker_from_to:bad_gain_dim', ...
              'Kp_tr and Kd_tr must be 2x2.');
    end

    % ---------- trailer velocity & error (desired vel = 0) ----------
    vxt = x(10);
    vyt = x(11);
    v_tr   = [vxt; vyt];
    v_tr_d = [0; 0];
    ev_tr  = v_tr - v_tr_d;

    % ---------- virtual planar force on trailer ----------
    f_xy = - Kp_tr * e_tr - Kd_tr * ev_tr;   % 2x1

    % ---------- map f_xy to [F_drive; tau_r] ----------
    theta_r = x(3);
    R = [cos(theta_r);  sin(theta_r)];   % heading direction
    T = [-sin(theta_r); cos(theta_r)];   % lateral direction

    F_nom   = dot(f_xy, R);             % along heading
    tau_nom = dot(f_xy, T);             % lateral → steering torque

    u_nom = [F_nom; tau_nom];

    % ---------- "QP": project onto box [u_min, u_max] ----------
    % 这是 min ||u - u_nom||^2 s.t. u_min <= u <= u_max 的解析解
    u = min(max(u_nom, u_min), u_max);

    if any(~isfinite(u))
        error('qp_path_tracker_from_to:NaN_in_u', ...
              'NaN/Inf in u at t=%.6f, u = [%g %g]', t, u(1), u(2));
    end

    % ---------- outputs ----------
    F_drive = u(1);
    tau_r   = u(2);

    % ---------- debug struct ----------
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
