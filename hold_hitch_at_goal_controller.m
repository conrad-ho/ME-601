function [F_drive, tau_r] = hold_hitch_at_goal_controller(x, t, params, y_goal)
% hold_hitch_at_goal_controller
%   在 hitch_ref 中用 'hold' 模式，把 hitch 固定在 y_goal 附近。
%   然后用 task_space_qp_controller 生成 [F_drive, tau_r].

    % 使用 persistent ref_state 让 hitch_ref 记住模式和 y0
    persistent ref_state

    if isempty(ref_state)
        ref_state = struct;
        ref_state.mode = 'hold';
        ref_state.y0   = y_goal;   % 把目标点设置成 hold 的位置
    end

    % 调用 hitch_ref 得到期望轨迹 (yd, ydot_d, yddot_ff)
    [yd, ydot_d, yddot_ff, ref_state] = hitch_ref(t, x, ref_state);

    % 再用你的 task-space QP 控制器求解
    [F_drive, tau_r, ~] = task_space_qp_controller(x, t, params, yd, ydot_d, yddot_ff);
end
