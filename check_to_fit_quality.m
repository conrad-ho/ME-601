function check_to_fit_quality(x_to, ref_to, params, dt_to)
%CHECK_TO_FIT_QUALITY 打印 TO 轨迹和参考路径的拟合误差
%
% x_to   : [nx x (N+1)] TO 计算得到的状态轨迹
%          state: [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
% ref_to : struct，包含字段
%           .t       [1 x (N+1)]
%           .p_hitch [2 x (N+1)]   hitch 参考路径
%           .p_tr    [2 x (N+1)]   trailer COM 参考路径
% params : 系统参数结构体（包含 d 等）
% dt_to  : 标称时间步长
    x_to=x_to.';
    % ---------- 基本尺寸检查 ----------
    [nx, Ncol] = size(x_to);
    N_ref = numel(ref_to.t);

    if Ncol ~= N_ref
        warning('x_to 列数 (%d) 与 ref_to.t 长度 (%d) 不一致！', Ncol, N_ref);
    end

    % 理论时间轴（如果你有真正的 tspan，也可以直接传进来用）
    t_theory = (0:Ncol-1) * dt_to;
    if isfield(ref_to, 't')
        dt_ref = diff(ref_to.t);
        fprintf('--- 时间步检查 ---\n');
        fprintf('  理论 dt_to      = %.4g\n', dt_to);
        fprintf('  ref_to dt mean  = %.4g, max dev = %.4g\n', ...
            mean(dt_ref), max(abs(dt_ref - dt_to)));
        fprintf('  t 对齐 max |t_theory - ref_to.t| = %.4g\n\n', ...
            max(abs(t_theory - ref_to.t)));
    end

    % ---------- 从 x_to 计算 hitch / trailer 位置 ----------
    p_hitch_to = zeros(2, Ncol);
    p_tr_to    = zeros(2, Ncol);

    for k = 1:Ncol
        xk = x_to(:,k);

        % Hitch：通过 helper 函数，从 robot pose + d 推出来
        p_hitch_to(:,k) = hitch_from_state(xk, params);  % 2x1

        % Trailer COM：state 里的 xt, yt
        % x = [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
        xt = xk(7);
        yt = xk(8);
        p_tr_to(:,k) = [xt; yt];
    end

    % ---------- Trailer 误差 ----------
    if isfield(ref_to, 'p_tr')
        e_tr = p_tr_to - ref_to.p_tr;
        e_tr_norm = vecnorm(e_tr);  % 每个时间点的欧氏距离误差

        max_err_tr = max(e_tr_norm);
        rms_err_tr = sqrt(mean(e_tr_norm.^2));

        fprintf('--- Trailer 轨迹拟合 ---\n');
        fprintf('  max |e_tr|  = %.4g (m)\n', max_err_tr);
        fprintf('  RMS |e_tr|  = %.4g (m)\n', rms_err_tr);
        
        % 起点终点误差
        fprintf('  start error = [%.4g, %.4g] (m)\n', e_tr(1,1), e_tr(2,1));
        fprintf('  end   error = [%.4g, %.4g] (m)\n\n', ...
                e_tr(1,end), e_tr(2,end));
    else
        warning('ref_to 中没有 p_tr 字段，无法计算 trailer 拟合误差。');
    end

    % ---------- Hitch 误差 ----------
    if isfield(ref_to, 'p_hitch')
        e_h = p_hitch_to - ref_to.p_hitch;
        e_h_norm = vecnorm(e_h);

        max_err_h = max(e_h_norm);
        rms_err_h = sqrt(mean(e_h_norm.^2));

        fprintf('--- Hitch 轨迹拟合 ---\n');
        fprintf('  max |e_h|   = %.4g (m)\n', max_err_h);
        fprintf('  RMS |e_h|   = %.4g (m)\n', rms_err_h);
        
        fprintf('  start error = [%.4g, %.4g] (m)\n', e_h(1,1), e_h(2,1));
        fprintf('  end   error = [%.4g, %.4g] (m)\n', ...
                e_h(1,end), e_h(2,end));
        fprintf('\n');
    else
        warning('ref_to 中没有 p_hitch 字段，无法计算 hitch 拟合误差。');
    end

    fprintf('==== TrajOpt vs ref_to 拟合检查完成 ====\n\n');
end
function p_tr = hitch_from_state(x,params) 
d = params.d; xr = x(1); yr = x(2); thetar = x(3); 
% hitch 点在世界坐标下的位置 
 p_tr = [xr - d*cos(thetar); 
 yr - d*sin(thetar)];
end