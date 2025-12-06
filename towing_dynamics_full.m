function xdot = towing_dynamics_full(x, u, params)
    import casadi.*

    % 拆 q,v
    q = [x(1:3); x(7:9)];
    v = [x(4:6); x(10:12)];

    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;      % 6x6
    B    = dyn.B;      % 6x2
    J    = dyn.J;      % 3x6
    dotJ = dyn.dotJ;   % 3x6
    vgen = dyn.v;      % 6x1

    % KKT 系统
    K   = [M, -J';
           J, SX.zeros(size(J,1), size(J,1))];    % 9x9
    rhs = [B*u;
          -dotJ*vgen];                            % 9x1

    sol = K \ rhs;
    a   = sol(1:6);   % generalized acceleration

    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end
