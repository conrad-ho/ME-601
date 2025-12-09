function debug_to_robot_trailer(x_to, ref_to, params, dt_to)
%DEBUG_TO_ROBOT_TRAILER
%   统一做几件事：
%   1) 把 x_to 解释成 robot+trailer 轨迹；
%   2) 和 ref_to.p_r / ref_to.p_tr 对齐，比 start/end error、max/RMS error；
%   3) 检查刚体几何：|p_tr - p_r| 是否接近 L_rt；
%   4) 画出 XY 轨迹 + 误差 + 几何误差。

    % -------- 形状处理：保证 x_to 是 [nx x (N+1)] --------
   x_to = x_to.';                % 现在 x_to: [nx x (N+1)]
    [nx, Ncol] = size(x_to); %#ok<NASGU>

    % 时间轴
    if isfield(ref_to,'t') && numel(ref_to.t)==Ncol
        t = ref_to.t(:).';
    else
        t = (0:Ncol-1)*dt_to;
    end

    % robot-trailer 几何参数
    if isfield(params,'d') && isfield(params,'Lt')
        L_rt = params.d + params.Lt;
    elseif isfield(params,'L_rt')
        L_rt = params.L_rt;
    else
        L_rt = NaN; % 暂时不知道
    end

    % -------- 从 x_to 提取 robot / trailer 轨迹 --------
    p_r_to  = zeros(2,Ncol);
    p_tr_to = zeros(2,Ncol);

    for k = 1:Ncol
        xk = x_to(:,k);
        % x = [xr, yr, thetar, vxr, vyr, wr,  xt, yt, thetat, vxt, vyt, wt]
        xr = xk(1);  yr = xk(2);
        xt = xk(7);  yt = xk(8);
        p_r_to(:,k)  = [xr; yr];
        p_tr_to(:,k) = [xt; yt];
    end

    % -------- 误差计算（注意这里不转置 ref_to） --------
    e_r  = []; e_tr = [];
    if isfield(ref_to,'p_r')
        e_r = p_r_to - ref_to.p_r;       % 都是 2 x N
    end
    if isfield(ref_to,'p_tr')
        e_tr = p_tr_to - ref_to.p_tr;    % 都是 2 x N
    end

    fprintf('=== TO vs ref (robot & trailer) ===\n');

    if ~isempty(e_tr)
        e_tr_norm = vecnorm(e_tr);
        fprintf('Trailer: max|e| = %.3f m, RMS|e| = %.3f m\n', ...
            max(e_tr_norm), sqrt(mean(e_tr_norm.^2)));
        fprintf('Trailer start error = [%.3f, %.3f] m\n', e_tr(1,1), e_tr(2,1));
        fprintf('Trailer end   error = [%.3f, %.3f] m\n', ...
            e_tr(1,end), e_tr(2,end));
    else
        fprintf('Trailer ref_to.p_tr 不存在，无法比较\n');
    end

    if ~isempty(e_r)
        e_r_norm = vecnorm(e_r);
        fprintf('Robot  : max|e| = %.3f m, RMS|e| = %.3f m\n', ...
            max(e_r_norm), sqrt(mean(e_r_norm.^2)));
        fprintf('Robot   start error = [%.3f, %.3f] m\n', e_r(1,1), e_r(2,1));
        fprintf('Robot   end   error = [%.3f, %.3f] m\n', ...
            e_r(1,end), e_r(2,end));
    else
        fprintf('Robot ref_to.p_r 不存在，无法比较\n');
    end

    % -------- 刚体几何检查：|p_tr - p_r| 是否恒定 --------
    if ~isnan(L_rt)
        dist_rt = vecnorm(p_tr_to - p_r_to);
        fprintf('Rigid geom: L_rt = %.3f m, mean|dist-L_rt| = %.4f m, max|dist-L_rt| = %.4f m\n', ...
            L_rt, mean(abs(dist_rt - L_rt)), max(abs(dist_rt - L_rt)));
    else
        fprintf('Rigid geom: L_rt 未知（params 里没有 d+Lt），跳过几何误差检查\n');
    end

    % -------- 图像：轨迹 & 误差 --------
    figure;
    tiledlayout(2,2);

    % (1) XY 轨迹
    nexttile;
    hold on; grid on; axis equal;
    if isfield(ref_to,'p_tr'), plot(ref_to.p_tr(1,:), ref_to.p_tr(2,:),'k--','LineWidth',1.5); end
    if isfield(ref_to,'p_r'),  plot(ref_to.p_r(1,:),  ref_to.p_r(2,:),'k:','LineWidth',1.5);   end
    plot(p_tr_to(1,:), p_tr_to(2,:),'b-','LineWidth',1.5);
    plot(p_r_to(1,:),  p_r_to(2,:),'r-','LineWidth',1.5);
    xlabel('x (m)'); ylabel('y (m)');
    title('XY trajectories: ref vs TO');
    legend({'ref trailer','ref robot','TO trailer','TO robot'},'Location','best');

    % (2) 误差随时间
    nexttile;
    hold on; grid on;
    if ~isempty(e_tr), plot(t, vecnorm(e_tr), 'b-', 'LineWidth',1.5); end
    if ~isempty(e_r),  plot(t, vecnorm(e_r),  'r--','LineWidth',1.5); end
    xlabel('time (s)'); ylabel('position error (m)');
    legend({'|e_{trailer}|','|e_{robot}|'},'Location','best');
    title('Tracking error vs time');

    % (3) 刚体距离
    if ~isnan(L_rt)
        nexttile;
        hold on; grid on;
        dist_rt = vecnorm(p_tr_to - p_r_to);
        plot(t, dist_rt, 'k-', 'LineWidth',1.2);
        yline(L_rt,'r--','LineWidth',1.2);
        xlabel('time (s)'); ylabel('|p_{tr} - p_r| (m)');
        title('Robot–trailer distance');
        legend({'actual','L\_rt target'},'Location','best');
    end

    % (4) 留空以后加别的 debug（比如 yaw）
end
