function ref_to = make_local_circle_ref(tspan, y0, opts)
% make_local_circle_ref
%   生成一个和 hitch_ref 里 'local_circle' 模式几乎一致的 hitch 参考轨迹，
%   用于给 towing_trajopt 的 ref_to.p_hitch 使用。
%
% inputs:
%   tspan - 1 x (N+1) 或 (N+1) x 1 的时间网格, 例如 0:0.02:30
%   y0    - 2x1, hitch 的初始位置 (对应 online 版里的 ref_state.y0)
%   opts  - 可选 struct:
%           opts.R         : 半径, 默认 1.0
%           opts.w         : 角速度, 默认 0.275
%           opts.phi0      : 初始角, 默认 -pi/2 (起始切线沿 +x)
%           opts.angle_max : 最大扫角, 默认 2*pi (绕一整圈)
%
% output:
%   ref_to.t        - 1 x (N+1) 时间网格 (复制自 tspan)
%   ref_to.p_hitch  - 2 x (N+1) hitch 参考位置 [xh; yh]
%   (兼容旧代码: ref_to.p_tr = p_hitch)

    if nargin < 3
        opts = struct;
    end
    if nargin < 2 || isempty(y0)
        y0 = [0; 0];
    end

    % ---- 确保 tspan 是行向量 ----
    tspan = tspan(:).';        % 强制成 1 x (N+1)
    Nt = numel(tspan);

    % ---- 读取 / 设定参数 (和 hitch_ref 保持一致) ----
    if ~isfield(opts, 'R'),         opts.R         = 1.0;    end
    if ~isfield(opts, 'w'),         opts.w         = 0.275;  end
    if ~isfield(opts, 'phi0'),      opts.phi0      = -pi/2;  end
    if ~isfield(opts, 'angle_max'), opts.angle_max = 2*pi;   end

    R         = opts.R;
    w         = opts.w;
    phi0      = opts.phi0;
    angle_max = opts.angle_max;

    % 圆心: 在 y0 上方 R
    center = y0 + [0; R];

    % 轨迹起始时间
    t0 = tspan(1);

    % 整个圆弧运动时间
    T_arc = angle_max / abs(w);

    % ---- 为每个时间点计算 theta(t) ----
    theta = zeros(1, Nt);
    for k = 1:Nt
        t = tspan(k);
        t_rel = t - t0;   % 相对起始时间

        if t_rel <= 0
            % 未开始, 停在起点
            theta_k = phi0;
        elseif t_rel <= T_arc
            % 圆弧运动中
            theta_k = phi0 + w * t_rel;
        else
            % 超出圆弧时间, 停在终点
            theta_k = phi0 + sign(w) * angle_max;
        end

        theta(k) = theta_k;
    end

    % ---- 根据 theta 计算 hitch 位置 ----
    % y = center + R * [cos(theta); sin(theta)]
    xh = center(1) + R * cos(theta);
    yh = center(2) + R * sin(theta);

    p_hitch = [xh; yh];   % 2 x (N+1)

    % ---- 输出 ref_to 结构 ----
    ref_to = struct;
    ref_to.t        = tspan;      % 1 x (N+1)
    ref_to.p_hitch  = p_hitch;    % 2 x (N+1)

    % 兼容旧代码: 如果 towing_trajopt 里还在用 ref.p_tr，可以顺手给一份
    ref_to.p_tr     = p_hitch;
end
