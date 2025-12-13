function dyn = towing_dynamics_mats(x, params)
    % x: 12x1 [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
    error('This Func is not been use, u are using wrong func');
    xr = x(1); yr = x(2); thetar = x(3);
    vxr = x(4); vyr = x(5); wr = x(6);
    xt = x(7); yt = x(8); thetat = x(9);
    vxt = x(10); vyt = x(11); wt = x(12);

    v = [vxr; vyr; wr; vxt; vyt; wt];

    m_r = params.m_r; I_r = params.I_r;
    m_t = params.m_t; I_t = params.I_t;
    d   = params.d;   Lt  = params.Lt;

    % 1) Mass Matrix
    M = diag([m_r m_r I_r m_t m_t I_t]);

    % 2) B
    B = [ cos(thetar), 0;
          sin(thetar), 0;
          0,           1;
          0,           0;
          0,           0;
          0,           0 ];

    % 3) constraint
    r_rh_body = [-d; 0];
    r_th_body = [Lt/2; 0];

    Rr = [cos(thetar) -sin(thetar);
          sin(thetar)  cos(thetar)];
    Rt = [cos(thetat) -sin(thetat);
          sin(thetat)  cos(thetat)];

    p_r = Rr * r_rh_body;
    p_t = Rt * r_th_body;

    S_p_r = [-p_r(2); p_r(1)];
    S_p_t = [-p_t(2); p_t(1)];

    % 4) constraint Jacobian J
    J_holo = [ eye(2), S_p_r, -eye(2), -S_p_t ];
    l = [-sin(thetat); cos(thetat)];
    J_nonholo = [0 0 0 l(1) l(2) 0];
    J = [J_holo; J_nonholo];

    % 5) dotJ
    Sd_p_r = [ d*cos(thetar)*wr; d*sin(thetar)*wr ];
    Sd_p_t = [ -(Lt/2)*cos(thetat)*wt; -(Lt/2)*sin(thetat)*wt ];

    dotJ_holo = [ zeros(2,2), Sd_p_r, zeros(2,2), -Sd_p_t ];
    ldot = [-cos(thetat)*wt; -sin(thetat)*wt];
    dotJ_nonholo = [0 0 0 ldot(1) ldot(2) 0];
    dotJ = [dotJ_holo; dotJ_nonholo];

    % 6) J and dJ
    J_y = [ 1, 0,  d*sin(thetar),  0, 0, 0;
            0, 1, -d*cos(thetar),  0, 0, 0 ];
    Jdot_y_v = [ d*cos(thetar)*wr;
                 d*sin(thetar)*wr ];

    % wrapper to struct
    dyn.M = M;
    dyn.B = B;
    dyn.J = J;
    dyn.dotJ = dotJ;
    dyn.v = v;
    dyn.J_y = J_y;
    dyn.Jdot_y_v = Jdot_y_v;
    dyn.H=zeros(6,1);
end
