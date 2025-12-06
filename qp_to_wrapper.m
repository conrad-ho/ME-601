function [F_drive, tau_r] = qp_to_wrapper(x, t, params)
% QP wrapper: 利用 TO 生成的 hitch 轨迹作为参考
% 要求 params.to_traj 已经由 build_to_qp_traj 构造好

    traj = params.to_traj;
    d    = params.d;

    % ---------- 当前 hitch 输出（仅用于 debug，可选） ----------
    dyn       = towing_dynamics_mats(x, params);
    J_y       = dyn.J_y;
    vgen      = dyn.v;

    xr     = x(1);
    yr     = x(2);
    thetar = x(3);
    y      = [xr - d*cos(thetar);
              yr - d*sin(thetar)];  %#ok<NASGU>

    % ---------- 从 traj_to 中插值出参考 yd, ydot_d, yddot_ff ----------
    t_ref     = traj.t;
    t_clamped = min(max(t, t_ref(1)), t_ref(end));   % 超出范围就 clamp

    yd        = interp1(t_ref, traj.y.'       , t_clamped, 'pchip').';
    ydot_d    = interp1(t_ref, traj.ydot.'    , t_clamped, 'pchip').';
    yddot_ff  = interp1(t_ref, traj.yddot_ff.', t_clamped, 'pchip').';

    % ---------- 调用原来的 task-space QP 控制器 ----------
    [F_drive, tau_r] = task_space_qp_controller( ...
                        x, t, params, yd, ydot_d, yddot_ff);
end
