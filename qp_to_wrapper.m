function [F_drive, tau_r] = qp_to_wrapper(x, t, params)
    traj = params.to_traj;
    d    = params.d;

    % 当前 hitch 输出（如果 task_space_qp_controller_proj 里不需要可以删掉）
    dyn  = towing_dynamics_mats(x, params);
    J_y  = dyn.J_y;
    vgen = dyn.v;

    xr     = x(1);
    yr     = x(2);
    thetar = x(3);
    y      = [xr - d*cos(thetar);
              yr - d*sin(thetar)];
    ydot   = J_y * vgen;  %#ok<NASGU>

    % ==== 用索引而不是 interp1 对齐 TO 轨迹 ====
    t_ref = traj.t;                  % 1 x (N+1)
    dt    = t_ref(2) - t_ref(1);     % assume uniform

    if t <= t_ref(1)
        k = 1;
    elseif t >= t_ref(end)
        k = numel(t_ref);
    else
        k = floor((t - t_ref(1))/dt) + 1;
        k = max(1, min(k, numel(t_ref)));
    end

    yd       = traj.y(:,       k);   % 2x1
    ydot_d   = traj.ydot(:,    k);   % 2x1
    yddot_ff = traj.yddot_ff(:,k);   % 2x1

    % ==== 调回原来的 CasADi 版 QP ====
    [F_drive, tau_r] = task_space_qp_controller_proj( ...
                        x, t, params, yd, ydot_d, yddot_ff);
end
