function [xdot, lambda] = towing_dynamics_full_with_lambda(x, u, params)
% towing_dynamics_full_with_lambda
%   基于 towing_dynamics_mats，解 KKT 得到
%   - generalized acceleration a
%   - 约束乘子 lambda（holo + nonholo）
%
% x : 12x1 state
% u :  2x1 control
% params : struct, 包含 m_r, m_t, I_r, I_t, d, Lt, 等

    import casadi.*

    % generalized coordinates / velocities（这里主要是为了 dyn.v）
    q = [x(1:3); x(7:9)];
    v = [x(4:6); x(10:12)]; %#ok<NASGU>

    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;      % 6x6
    B    = dyn.B;      % 6x2
    J    = dyn.J;      % 3x6
    dotJ = dyn.dotJ;   % 3x6
    vgen = dyn.v;      % 6x1

    % KKT 系统
    K   = [M, -J';
           J, SX.zeros(size(J,1), size(J,1))];
    rhs = [B*u;
          -dotJ*vgen];

    sol = K \ rhs;
    a       = sol(1:6);                 % generalized acceleration
    lambda  = sol(7:end);               % 3x1 constraint forces

    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end
%% ======== helper functions (保持你现有的定义) ========

function xdot = towing_dynamics_full(x, u, params)
    import casadi.*

    q = [x(1:3); x(7:9)];
    v = [x(4:6); x(10:12)];

    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;
    B    = dyn.B;
    J    = dyn.J;
    dotJ = dyn.dotJ;
    vgen = dyn.v;

    K   = [M, -J';
           J, SX.zeros(size(J,1), size(J,1))];
    rhs = [B*u;
          -dotJ*vgen];

    sol = K \ rhs;
    a   = sol(1:6);

    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end

function p_h = hitch_from_state(x,params)
    d  = params.d;
    xr = x(1);  yr = x(2);  thetar = x(3);
    p_h = [xr - d*cos(thetar);
           yr - d*sin(thetar)];
end

function phi = hitch_angle_from_state(x)
    theta_r = x(3);
    theta_t = x(9);
    phi = theta_t - theta_r;
end
