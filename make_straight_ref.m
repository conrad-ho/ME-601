function ref_to = make_straight_ref(tspan, y0, opts)
% make_S_arc_ref
%   生成一个简单的 S 形 hitch 参考路径:
%   先左转一个固定角度, 再右转同样角度, 半径相同.
%
% inputs:
%   tspan : 1x(N+1) 或 (N+1)x1 时间网格
%   y0    : 2x1 起始 hitch 位置
%   opts  : 选项
%           opts.R      : 圆弧半径, 默认 1.0
%           opts.w      : 角速度, 默认 0.3 (rad/s)
%           opts.angle  : 每段扫角, 默认 pi/2 (90 度)
%           opts.theta0 : 初始朝向, 默认 0 (沿 +x)
%
% output:
%   ref_to.t       : 1x(N+1) 时间
%   ref_to.p_hitch : 2x(N+1) [xh; yh]
%   ref_to.p_tr    : 同 p_hitch, 为兼容旧代码

    if nargin < 3
        opts = struct;
    end
    if nargin < 2 || isempty(y0)
        y0 = [0;0];
    end

    if ~isfield(opts,'R'),      opts.R      = 2.0;     end
    if ~isfield(opts,'w'),      opts.w      = 0.3;     end
    if ~isfield(opts,'angle'),  opts.angle  = pi/2;    end
    if ~isfield(opts,'theta0'), opts.theta0 = 0.0;     end  % 初始 heading

    R      = opts.R;
    w      = opts.w;              % 角速度 (左转为 +w)
    angle  = opts.angle;          % 每段扫角
    theta0 = opts.theta0;

    % 速度 v = R * |w|, 这样曲率刚好是 1/R
    v = R * abs(w);

    % 每段时间
    T1 = angle / abs(w);    % 左转时长
    T2 = 2*T1;              % 左+右结束时间

    % 处理时间网格
    tspan = tspan(:).';
    Nt    = numel(tspan);
    dt    = tspan(2) - tspan(1);  % 假设等间距

    % 初始化
    p     = zeros(2, Nt);
    theta = zeros(1, Nt);

    p(:,1)     = y0(:);     % 位置
    theta(1)   = theta0;    % 朝向

    for k = 2:Nt
        t_prev = tspan(k-1);
        th     = theta(k-1);

        % 当前时间落在哪一段
        if t_prev <= T1
            % 第一段: 左转 (曲率 +1/R)
            w_k = sign(w) * abs(w);
        elseif t_prev <= T2
            % 第二段: 右转 (曲率 -1/R)
            w_k = -sign(w) * abs(w);
        else
            % 后面保持直行 (或你也可以设 v=0 就原地停)
            w_k = 0;
        end

        % 积分朝向
        theta(k) = th + w_k * dt;

        % 积分位置: p_dot = v * [cos(theta); sin(theta)]
        p(:,k) = p(:,k-1) + v * [cos(th); sin(th)] * dt;
    end

    p_hitch = p;

    ref_to = struct;
    ref_to.t       = tspan;
    ref_to.p_hitch = p_hitch;
    ref_to.p_tr    = p_hitch;   % 兼容旧代码
end
