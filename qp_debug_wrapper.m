function [F_drive, tau_r] = qp_debug_wrapper(x, t, params)
% wrapper around task_space_qp_controller for debugging
% interface: [F_drive, tau_r] = qp_debug_wrapper(x, t, params)

    d = params.d;
    
    % -------- current hitch output y, ydot (for debugging and reference) --------
    % state unpacking: according to your system definition, assume x = [q(1:6); dq(1:6)]
    xr     = x(1);
    yr     = x(2);
    thetar = x(3);
   

    y = [xr - d*cos(thetar);
         yr - d*sin(thetar)];

    dyn       = towing_dynamics_mats(x, params);
    v      = dyn.v;      % change this if dq is not located in 7:12
    J_y       = dyn.J_y;
    Jdot_y_v  = dyn.Jdot_y_v;
    ydot      = J_y * v;

    persistent ref_state
    if isempty(ref_state)
        ref_state.mode = 'L_turn';   % <<<<<< change reference mode here if needed
        ref_state.initialized = false;
    end
    %ref_state.mode = 'hold';
    %ref_state.mode = 'step_x';
    ref_state.mode = 'local_circle';
    %ref_state.mode = 'L_turn';

    [yd, ydot_d, yddot_ff, ref_state] = hitch_ref(t, y, ref_state);

    % ======================================================
    % 2) call the QP controller: note we must pass yd, ydot_d, yddot_ff
    % ======================================================
    %[F_drive, tau_r, qp_dbg] = task_space_qp_controller( ...
    %                        x, t, params, yd, ydot_d, yddot_ff);
    [F_drive, tau_r, qp_dbg] = task_space_qp_controller( ...
                          x, t, params, yd, ydot_d, yddot_ff);
    % unpack qp solution
    a_opt      = qp_dbg.a_opt;
    u_opt      = qp_dbg.u_opt;
    lambda_opt = qp_dbg.lambda_opt;
    exitflag   = qp_dbg.exitflag;

    % ======================================================
    % 3) compute constraint residuals and tracking metrics (optional)
    % ======================================================

    % re-fetch dynamics matrices (dyn is already computed above)
    M    = dyn.M;
    B    = dyn.B;
    J    = dyn.J;
    dotJ = dyn.dotJ;
    v    = dyn.v;          % use v from dyn for consistency

    % equality constraint residuals
    dyn_res = M * a_opt - B * u_opt - J' * lambda_opt;   % should be close to 0
    con_res = J * a_opt + dotJ * v;                      % should be close to 0

    % task-space acceleration tracking error
    Kp = diag([25, 25]);     % must match the controller implementation
    Kd = diag([1.5, 1.5]);

    yddot_des = yddot_ff ...
                - Kd * (ydot - ydot_d) ...
                - Kp * (y    - yd);

    yddot_actual = J_y * a_opt + Jdot_y_v;
    yddot_err    = yddot_actual - yddot_des;

    % simple printout (comment out if too verbose)
    % fprintf('t = %.3f | exit = %2d | ||dyn_res|| = %.2e | ||con_res|| = %.2e | ||yddot_err|| = %.2e | u = [%.3f, %.3f]\n', ...
    %         t, exitflag, norm(dyn_res), norm(con_res), norm(yddot_err), u_opt(1), u_opt(2));

end
