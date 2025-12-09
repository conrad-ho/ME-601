function ref_to = make_to_ref_Sarc_robot_trailer(tspan, y0, opts)
% make_to_ref_Sarc_robot_trailer
%   为 TO 生成一条简单的 S 形参考路径，
%   输出 trailer 参考轨迹 p_tr、hitch 轨迹 p_h，以及与之满足刚体约束的 robot 参考轨迹 p_r。
%
% 几何约定：
%   - theta(k) 表示 “从 trailer 指向 robot” 的朝向（也是 robot 朝向）。
%   - p_tr(:,k)  : trailer 质心位置
%   - p_h(:,k)   : hitch 位置（trailer 前方）
%   - p_r(:,k)   : robot 质心位置（在 hitch 前方）
%
%   也就是说，在每个 k：
%       n_k = [cos(theta(k)); sin(theta(k))]   % 单位朝向向量（trailer -> robot）
%       p_h(:,k) = p_tr(:,k) + Lt * n_k
%       p_r(:,k) = p_h(:,k) + d  * n_k
%       L_rt     = Lt + d      = norm(p_r - p_tr) （名义距离）
%
% inputs:
%   tspan : 1x(N+1) 或 (N+1)x1 时间网格（等间距）
%   y0    : 2x1 trailer 初始位置 p_tr(:,1)
%   opts  : 选项结构体
%           opts.R       : 圆弧半径, 默认 2.0
%           opts.w       : 名义角速度, 默认 0.3 (rad/s)
%           opts.angle   : 每段扫角, 默认 pi/2 (90 度)
%           opts.theta0  : 初始朝向（trailer->robot）, 默认 0 (沿 +x)
%
%   几何参数（必须能推出 Lt 和 d）：
%     方案 1（推荐、和动力学一致）：
%           opts.d       : robot COM -> hitch 的距离（沿 theta 方向）
%           opts.Lt      : trailer COM -> hitch 的距离（沿 theta 方向）
%           ==> L_rt = d + Lt
%     方案 2：
%           opts.L_rt    : robot COM -> trailer COM 的距离
%           opts.Lt      : trailer COM -> hitch 的距离
%           ==> d = L_rt - Lt（若 d < 0 则会报错）
%
% output:
%   ref_to.t     : 1x(N+1) 时间
%   ref_to.p_tr  : 2x(N+1) trailer 参考位置 [x_t; y_t]
%   ref_to.p_h   : 2x(N+1) hitch   参考位置 [x_h; y_h]
%   ref_to.p_r   : 2x(N+1) robot   参考位置 [x_r; y_r]
%   ref_to.theta : 1x(N+1) 参考朝向（trailer -> robot）

    % --------- 默认参数处理 ---------
    if nargin < 3
        opts = struct;
    end
    if nargin < 2 || isempty(y0)
        y0 = [0;0];
    end

    if ~isfield(opts,'R'),      opts.R      = 2.0;   end
    if ~isfield(opts,'w'),      opts.w      = 0.3;   end
    if ~isfield(opts,'angle'),  opts.angle  = pi/2;  end
    if ~isfield(opts,'theta0'), opts.theta0 = 0.0;   end  % 初始 heading (trailer -> robot)

    % --------- 解析几何参数 (Lt, d, L_rt) ---------
    has_d   = isfield(opts,'d');
    has_Lt  = isfield(opts,'Lt');
    has_Lrt = isfield(opts,'L_rt');

    if has_d && has_Lt
        % 推荐方式：直接给 d 和 Lt
        d    = opts.d;
        Lt   = opts.Lt/2; %这里Lt=COM-> hitch
        L_rt = d + Lt;
    elseif has_Lrt && has_Lt
        % 次选方式：给 L_rt 和 Lt，自动算 d
        L_rt = opts.L_rt;
        Lt   = opts.Lt;
        d    = L_rt - Lt;
        if d < -1e-9
            error('make_to_ref_Sarc_robot_trailer: L_rt (%.4g) 必须大于等于 Lt (%.4g)，否则几何不合理。', ...
                  L_rt, Lt);
        end
    else
        error(['make_to_ref_Sarc_robot_trailer: 需要几何参数以定义 trailer-hitch-robot 关系。\n' ...
               '请提供 (opts.d 和 opts.Lt)，或者 (opts.L_rt 和 opts.Lt)。']);
    end

    R      = opts.R;
    w_nom  = opts.w;         % 名义角速度 (左转为 +w_nom)
    angle  = opts.angle;     % 每段扫角
    theta0 = opts.theta0;

    % 速度 v = R * |w|，这样曲率刚好是 1/R
    v = R * abs(w_nom);

    % 每段时间
    T1 = angle / abs(w_nom);    % 第一段左转时长
    T2 = 2*T1;                  % 左+右结束时间

    % --------- 时间网格 ---------
    tspan = tspan(:).';
    Nt    = numel(tspan);
    if Nt < 2
        error('tspan 至少需要包含两个时间点');
    end
    dt = tspan(2) - tspan(1);  % 假设等间距

    % --------- 初始化 ---------
    p_tr  = zeros(2, Nt);   % trailer 位置
    theta = zeros(1, Nt);   % 参考朝向（trailer -> robot）

    p_tr(:,1) = y0(:);      % trailer 初始位置
    theta(1)  = theta0;     % 初始 heading

    % --------- 生成 trailer 轨迹 (S 形) ---------
    for k = 2:Nt
        t_prev = tspan(k-1);
        th     = theta(k-1);

        % 当前时间落在哪一段
        if t_prev <= T1
            % 第一段: 左转 (曲率 +1/R)
            w_k = sign(w_nom) * abs(w_nom);
        elseif t_prev <= T2
            % 第二段: 右转 (曲率 -1/R)
            w_k = -sign(w_nom) * abs(w_nom);
        else
            % 后面保持直行
            w_k = 0;
        end

        % 积分朝向
        theta(k) = th + w_k * dt;

        % 积分 trailer 位置: p_tr_dot = v * [cos(theta); sin(theta)]
        p_tr(:,k) = p_tr(:,k-1) + v * [cos(th); sin(th)] * dt;
    end

    % --------- 根据几何关系生成 hitch 和 robot 轨迹 ---------
    % 向量化：n = [cos(theta); sin(theta)] (2 x Nt)
    n_x = cos(theta);
    n_y = sin(theta);
    n   = [n_x; n_y];

    % hitch: 在 trailer 前方 Lt
    p_h = p_tr + Lt  * n;

    % robot: 在 hitch 前方 d （或 equivalently 在 trailer 前方 L_rt）
    p_r = p_tr + L_rt * n;
    % 等价：p_r = p_h + d * n;

    % --------- 初始 / 结束 几何一致性检查 ---------
    tol_init  = 1e-8;   % 初始要求非常严格
    tol_final = 1e-6;   % 结束允许一点数值误差

    idx_list = [1, Nt];
    for idx = idx_list
        if idx == 1
            tol = tol_init;
            tag = 'initial';
        else
            tol = tol_final;
            tag = 'final';
        end

        th_k  = theta(idx);
        n_k   = [cos(th_k); sin(th_k)];
        p_trk = p_tr(:,idx);
        p_hk  = p_h(:,idx);
        p_rk  = p_r(:,idx);

        % 三个约束：tr->h, h->r, tr->r
        e_tr_h = norm(p_hk - (p_trk + Lt  * n_k));
        e_h_r  = norm(p_rk - (p_hk  + d   * n_k));
        e_tr_r = norm(p_rk - (p_trk + L_rt * n_k));

        e_max = max([e_tr_h, e_h_r, e_tr_r]);

        if e_max > tol || any(~isfinite([e_tr_h, e_h_r, e_tr_r]))
            error(['make_to_ref_Sarc_robot_trailer: %s geometry check failed (err = %.3g).\n' ...
                   '请检查 tspan / theta / 几何参数 (d, Lt, L_rt) 是否与动力学约定一致。'], ...
                   tag, e_max);
        end
    end

    % --------- 打包输出 ---------
    ref_to = struct;
    ref_to.t     = tspan;   % 1x(N+1)
    ref_to.p_tr  = p_tr;    % 2x(N+1)
    ref_to.p_h   = p_h;     % 2x(N+1)
    ref_to.p_r   = p_r;     % 2x(N+1)
    ref_to.theta = theta;   % 1x(N+1)
end
