function traj_to = build_to_qp_traj(dt, x_sol, u_sol, params)
% build_to_qp_traj
% build hitch reference for QP from TO solution (x_sol, u_sol):
%   traj_to.t, traj_to.y, traj_to.ydot, traj_to.yddot_ff
%
% now uses unified dynamics interface:
%   towing_dynamic(x, params, need)

    N = size(u_sol,1);  % number of control steps

    traj_to = struct;
    traj_to.t = (0:N) * dt;   % 1 x (N+1)
    traj_to.x = x_sol;        % (N+1) x 12
    traj_to.u = u_sol;        % N x 2

    traj_to.y        = zeros(2, N+1);
    traj_to.ydot     = zeros(2, N+1);
    traj_to.yddot_ff = zeros(2, N+1);

    % compute y and ydot at all nodes
    for k = 1:N+1
        xk = x_sol(k,:).';

        dyn_m = towing_dynamic(xk, params, 'mats');
        traj_to.y(:,k)    = dyn_m.y;      % hitch position (2x1)
        traj_to.ydot(:,k) = dyn_m.ydot;   % hitch velocity  (2x1)
    end

    % compute feedforward yddot at nodes 1..N using task-space affine dynamics
    for k = 1:N
        xk = x_sol(k,:).';
        uk = u_sol(k,:).';

        dyn_task = towing_dynamic(xk, params, 'task');
        A_y = dyn_task.A_y;   % 2x2
        b_y = dyn_task.b_y;   % 2x1

        traj_to.yddot_ff(:,k) = A_y * uk + b_y;
    end

    % last yddot_ff: repeat final available
    traj_to.yddot_ff(:,N+1) = traj_to.yddot_ff(:,N);
end
