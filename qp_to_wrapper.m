function [F_drive, tau_r, qp_dbg] = qp_to_wrapper(x, t, params)
% qp_to_wrapper
% QP wrapper: uses TO-generated hitch trajectory as reference
% requires params.to_traj built by build_to_qp_traj

    traj = params.to_traj;

    % ---------- current hitch output (optional debug) ----------
    dyn_m = towing_dynamic(x, params, 'mats');
    % y    = dyn_m.y;      %#ok<NASGU>
    % ydot = dyn_m.ydot;   %#ok<NASGU>
    % J_y  = dyn_m.J_y;    %#ok<NASGU>
    % vgen = dyn_m.v;      %#ok<NASGU>

    % ---------- interpolate reference yd, ydot_d, yddot_ff from traj ----------
    t_ref     = traj.t(:);  % (N+1)x1
    t_clamped = min(max(t, t_ref(1)), t_ref(end));

    yd       = interp1(t_ref, traj.y.'       , t_clamped, 'pchip').';
    ydot_d   = interp1(t_ref, traj.ydot.'    , t_clamped, 'pchip').';
    yddot_ff = interp1(t_ref, traj.yddot_ff.', t_clamped, 'pchip').';

    % ---------- call projected task-space QP controller ----------
    [F_drive, tau_r, qp_dbg] = task_space_qp_controller_proj( ...
                        x, t, params, yd, ydot_d, yddot_ff);
end
