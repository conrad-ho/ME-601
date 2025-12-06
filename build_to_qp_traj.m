function traj_to = build_to_qp_traj(dt, x_sol, u_sol, params)
% 根据 TO 的离散解 (x_sol, u_sol)，构造给 QP 用的 hitch 参考：
%   traj_to.t, traj_to.y, traj_to.ydot, traj_to.yddot_ff

    N  = size(u_sol,1);   % 控制步数

    traj_to = struct;
    traj_to.t = (0:N) * dt;   % 1 x (N+1)
    traj_to.x = x_sol;        % (N+1) x 12
    traj_to.u = u_sol;        % N x 2

    traj_to.y        = zeros(2, N+1);
    traj_to.ydot     = zeros(2, N+1);
    traj_to.yddot_ff = zeros(2, N+1);

    d = params.d;

    % 先计算所有节点的 y, ydot
    for k = 1:N+1
        xk = x_sol(k,:).';

        % hitch 位置
        xr     = xk(1);
        yr     = xk(2);
        thetar = xk(3);
        traj_to.y(:,k) = [xr - d*cos(thetar);
                          yr - d*sin(thetar)];

        dyn  = towing_dynamics_mats(xk, params);
        J_y  = dyn.J_y;
        vgen = dyn.v;          % generalized velocity
        traj_to.ydot(:,k) = J_y * vgen;
    end

    % 再计算 0..N-1 节点处的 feedforward yddot
    for k = 1:N
        xk = x_sol(k,:).';
        uk = u_sol(k,:).';

        dyn  = towing_dynamics_mats(xk, params);
        M    = dyn.M;          % 6x6
        B    = dyn.B;          % 6x2
        J    = dyn.J;          % 3x6
        dotJ = dyn.dotJ;       % 3x6
        vgen = dyn.v;          % 6x1
        J_y  = dyn.J_y;        % 2x6
        Jdot_y_v = dyn.Jdot_y_v; % 2x1

        % KKT 系统解 generalized acceleration a_k
        K   = [M, -J';
               J, zeros(size(J,1))];      % 9x9
        rhs = [B*uk;
              -dotJ*vgen];                % 9x1

        sol = K \ rhs;
        a   = sol(1:6);                   % generalized acceleration

        % hitch 加速度
        yddot = J_y * a + Jdot_y_v;       % 2x1
        traj_to.yddot_ff(:,k) = yddot;
    end

    % 末端的 yddot_ff 就直接重复最后一个
    traj_to.yddot_ff(:,N+1) = traj_to.yddot_ff(:,N);
end
