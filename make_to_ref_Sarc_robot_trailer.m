function ref_to = make_to_ref_Sarc_robot_trailer(tspan, y0, opts)
% make_to_ref_Sarc_robot_trailer
%   为 TO 生成一条简单的 S 形参考路径，
%   输出 trailer 参考轨迹 p_tr 以及与之满足刚体约束的 robot 参考轨迹 p_r。
%
%   函数本质上还是沿用原来的“左转再右转”的 S 形轨迹，只是：
%   - y0 现在表示 trailer 的初始位置；
%   - p_tr 是 trailer 的路径；
%   - p_r 通过 kinematic 关系从 p_tr 和 heading 计算得到。
%
% inputs:
%   tspan : 1x(N+1) 或 (N+1)x1 时间网格（等间距）
%   y0    : 2x1 trailer 初始位置 p_tr(:,1)
%   opts  : 选项结构体
%           opts.R       : 圆弧半径, 默认 2.0
%           opts.w       : 角速度, 默认 0.3 (rad/s)
%           opts.angle   : 每段扫角, 默认 pi/2 (90 度)
%           opts.theta0  : 初始朝向, 默认 0 (沿 +x)
%
%           下面是 robot-trailer 几何参数（选其中一种方式给）：
%           方式 A：直接给 robot 到 trailer COM 的距离
%           opts.L_rt    : robot COM -> trailer COM 的距离
%
%           方式 B：如果你已经有 d 和 Lt（和动力学一致）：
%           opts.d       : robot COM -> hitch 的距离
%           opts.Lt      : hitch -> trailer COM 的距离
%           会自动用 L_rt = d + Lt
%
% output:
%   ref_to.t    : 1x(N+1) 时间
%   ref_to.p_tr : 2x(N+1) trailer 参考位置 [x_t; y_t]
%   ref_to.p_r  : 2x(N+1) robot   参考位置 [x_r; y_r]
%   ref_to.theta: 1x(N+1) 参考朝向（这里假设 robot 和 trailer 朝向相同）

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
    if ~isfield(opts,'theta0'), opts.theta0 = 0.0;   end  % 初始 heading

    % 处理 robot-trailer 距离 L_rt
    if isfield(opts,'d') && isfield(opts,'Lt')
        L_rt = opts.d + opts.Lt;     % 和动力学保持一致
    elseif isfield(opts,'L_rt')
        L_rt = opts.L_rt;
    else
        L_rt = 0.5;                  % 没给的话用一个温和的默认值
    end

    R      = opts.R;
    w_nom  = opts.w;                 % 名义角速度 (左转为 +w_nom)
    angle  = opts.angle;             % 每段扫角
    theta0 = opts.theta0;

    % 速度 v = R * |w|，这样曲率刚好是 1/R
    v = R * abs(w_nom);

    % 每段时间
    T1 = angle / abs(w_nom);    % 左转时长
    T2 = 2*T1;                  % 左+右结束时间

    % --------- 时间网格 ---------
    tspan = tspan(:).';
    Nt    = numel(tspan);
    if Nt < 2
        error('tspan 至少需要包含两个时间点');
    end
    dt    = tspan(2) - tspan(1);  % 假设等间距

    % --------- 初始化 ---------
    p_tr   = zeros(2, Nt);   % trailer 位置
    theta  = zeros(1, Nt);   % 参考朝向（这里假设 robot 和 trailer 相同）

    p_tr(:,1)  = y0(:);      % trailer 初始位置
    theta(1)   = theta0;     % 初始 heading

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

    % --------- 根据几何关系生成 robot 轨迹 ---------
    % 假设 robot 与 trailer 共线，并在前方（或后方）距离 L_rt
    % p_r = p_tr - L_rt * [cos(theta); sin(theta)]
    p_r = zeros(2,Nt);
    for k = 1:Nt
        th = theta(k);
        p_r(:,k) = p_tr(:,k) - L_rt * [cos(th); sin(th)];
    end

    % --------- 打包输出 ---------
    ref_to = struct;
    ref_to.t     = tspan;   % 1x(N+1)
    ref_to.p_tr  = p_tr;    % 2x(N+1)
    ref_to.p_r   = p_r;     % 2x(N+1)
    ref_to.theta = theta;   % 1x(N+1)，方便 TO 里加 yaw 相关 cost/constraint
end
