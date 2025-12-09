function check_kinematic_constraint_x_to(x_to, params, dt)
% check_kinematic_constraint_x_to
%   在 TO 结果 x_to 上检查 robot-trailer 刚体约束是否在几何上被保持。
%
% inputs:
%   x_to   : (N+1) x 12, towing_trajopt 输出的状态轨迹
%   params : 含有 d, Lt 等几何参数的 struct
%   dt     : 时间步长（仅用于画图，如果不想画图可以省略）

    if nargin < 3
        dt = 1.0; % 没给就用 1，当成归一化时间
    end

    [Nplus1, nx] = size(x_to);
    if nx ~= 12
        error('x_to 维度不对，期望 (N+1) x 12，得到 (N+1) x %d', nx);
    end

    N = Nplus1 - 1;
    t = (0:N) * dt;

    % 预分配误差
    e_pos   = zeros(Nplus1,1);   % 几何约束：两种算法的 hitch 位置差
    e_Jv    = zeros(Nplus1,1);   % 速度层约束：J*v 是否为 0（可选，需要 towing_dynamics_mats）

    have_dyn = exist('towing_dynamics_mats','file') == 2;

    for k = 1:Nplus1
        xk = x_to(k,:).';  % 12x1

        % --- 1) 从 robot 状态计算 hitch 位置 ---
        xr     = xk(1);
        yr     = xk(2);
        thetar = xk(3);
        d      = params.d;       % robot COM 到 hitch 的距离

        p_h_from_robot = [xr - d*cos(thetar);
                          yr - d*sin(thetar)];

        % --- 2) 从 trailer 状态计算 hitch 位置 ---
        xt      = xk(7);
        yt      = xk(8);
        thetat  = xk(9);
        Lt      = params.Lt;     % trailer 车体长度（假设 hitch 在前端）

        % 这里假设 hitch 在 trailer 前端，沿 heading 正方向 Lt/2
        % 如果你定义的是在后端，就改成 xt - (Lt/2)*cos(thetat) 等
        p_tr_com   = [xt; yt];
        p_h_from_tr = p_tr_com + (Lt/2) * [cos(thetat);
                                           sin(thetat)];

        % --- 3) 几何误差：两种计算方式的 hitch 位置差 ---
        e_pos(k) = norm(p_h_from_robot - p_h_from_tr);

        % --- 4) （可选）速度层约束：J*v 是否为 0 ---
        if have_dyn
            dyn = towing_dynamics_mats(xk, params);
            J   = dyn.J;     % 3 x 6
            v   = dyn.v;     % 6 x 1 (generalized速度)
            e_Jv(k) = norm(J * v);
        else
            e_Jv(k) = NaN;
        end
    end

    % ----- 打印统计量 -----
    fprintf('=== Kinematic constraint check (position level) ===\n');
    fprintf('  max ||p_h^robot - p_h^trailer|| = %.3e (m)\n', max(e_pos));
    fprintf('  RMS ||p_h^robot - p_h^trailer|| = %.3e (m)\n', rms(e_pos));

    if have_dyn
        fprintf('=== Velocity-level constraint check (J*v ≈ 0) ===\n');
        fprintf('  max ||J*v|| = %.3e\n', max(e_Jv));
        fprintf('  RMS ||J*v|| = %.3e\n', rms(e_Jv));
    else
        fprintf('未找到 towing_dynamics_mats.m，跳过 J*v 约束检查。\n');
    end

    % ----- 简单画图 -----
    figure;
    subplot(2,1,1);
    plot(t, e_pos, 'LineWidth', 1.5);
    grid on;
    xlabel('t');
    ylabel('||\Delta p_h|| (m)');
    title('几何层：robot vs trailer 计算的 hitch 位置误差');

    subplot(2,1,2);
    if have_dyn
        plot(t, e_Jv, 'LineWidth', 1.5);
        ylabel('||J v||');
        title('速度层：约束残差 ||J v||');
    else
        plot(t, e_pos*NaN);
        ylabel('N/A');
        title('未检查 Jv，因为缺少 towing_dynamics_mats');
    end
    xlabel('t');

end
    