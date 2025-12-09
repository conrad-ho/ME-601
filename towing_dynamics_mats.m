function dyn = towing_dynamics_mats(x, params)
%TOWING_DYNAMICS_MATS
%   Common robot+trailer matrices used by TO, simulation and QP.
%
% State:
%   x = [xr; yr; thetar; vxr; vyr; wr;
%        xt; yt; thetat; vxt; vyt; wt]  (12x1)
%
% Params must contain:
%   m_r, I_r, m_t, I_t, d, Lt

    % unpack state
    xr = x(1);  yr = x(2); thetar  = x(3);
    vxr = x(4); vyr = x(5); wr     = x(6);
    xt = x(7);  yt = x(8); thetat  = x(9);
    vxt = x(10); vyt = x(11); wt   = x(12);

    % generalized velocity
    v = [vxr; vyr; wr; vxt; vyt; wt];

    % params
    m_r = params.m_r; I_r = params.I_r;
    m_t = params.m_t; I_t = params.I_t;
    d   = params.d;   Lt  = params.Lt;

    % === 1) Mass matrix (6x6) ===
    M = diag([m_r m_r I_r m_t m_t I_t]);

    % === 2) Input matrix B (6x2) ===
    % u = [F_drive; tau_r], F_drive is along robot body +x
    cr = cos(thetar); sr = sin(thetar);
    B = [ cr, 0;
          sr, 0;
          0,  1;
          0,  0;
          0,  0;
          0,  0 ];

    % === 3) Constraint Jacobian J (3x6) ===
    % hitch offsets in body frames
    r_rh_body = [-d;   0];   % robot COM -> hitch
    r_th_body = [ Lt/2; 0];  % trailer COM -> hitch

    % rotation matrices
    ct = cos(thetat); st = sin(thetat);
    Rr = [ cr, -sr;
           sr,  cr];
    Rt = [ ct, -st;
           st,  ct];

    % hitch vectors in world frame
    p_r = Rr * r_rh_body;   % 2x1
    p_t = Rt * r_th_body;   % 2x1

    % skew-like terms for angular velocity contribution
    S_p_r = [-p_r(2); p_r(1)];
    S_p_t = [-p_t(2); p_t(1)];

    % holonomic: hitch points coincide
    % [vx_r; vy_r] + S_p_r*wr  - [vx_t; vy_t] - S_p_t*wt = 0
    J_holo = [ eye(2),  S_p_r,  -eye(2), -S_p_t ];   % 2x6

    % nonholonomic: trailer lateral velocity = 0 (body y-axis)
    l = [-sin(thetat); cos(thetat)];                 % body y in world
    J_nonholo = [0 0 0 l(1) l(2) 0];                 % 1x6

    J = [J_holo;
         J_nonholo];                                 % 3x6

    % === 4) Constraint Jacobian time derivative dotJ (3x6) ===
    % NOTE: keep same form as旧代码，避免引入新的几何差异
    Sd_p_r = [ d*cos(thetar)*wr;
               d*sin(thetar)*wr ];

    Sd_p_t = [ -(Lt/2)*cos(thetat)*wt;
                (Lt/2)*sin(thetat)*wt ];

    dotJ_holo = [ zeros(2,2), Sd_p_r, zeros(2,2), -Sd_p_t ];

    ldot = [-cos(thetat)*wt;
            -sin(thetat)*wt];
    dotJ_nonholo = [0 0 0 ldot(1) ldot(2) 0];

    dotJ = [dotJ_holo;
            dotJ_nonholo];

    % === 5) Hitch task Jacobian J_y (2x6) and Jdot_y_v (2x1) ===
    J_y = [ 1, 0,  d*sin(thetar),  0, 0, 0;
            0, 1, -d*cos(thetar),  0, 0, 0 ];

    Jdot_y_v = [ d*cos(thetar)*wr;
                 d*sin(thetar)*wr ];

    % pack struct
    dyn.M = M;
    dyn.B = B;
    dyn.J = J;
    dyn.dotJ = dotJ;
    dyn.v = v;
    dyn.J_y = J_y;
    dyn.Jdot_y_v = Jdot_y_v;
    dyn.H = zeros(6,1);  % placeholder if needed
end
