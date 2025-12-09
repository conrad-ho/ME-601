function xdot = towing_dynamics_full(x, u, params)
import casadi.*
%TOWING_DYNAMICS_FULL
%   Continuous-time rigid robot+trailer dynamics with Baumgarte-stabilized
%   holonomic (hitch) + nonholonomic constraints.
%
% State:
%   x = [xr; yr; thetar; vxr; vyr; wr;
%        xt; yt; thetat; vxt; vyt; wt]  (12x1)
%
% Input:
%   u = [F_drive; tau_r]
%
% Output:
%   xdot, same ordering as x: [vxr; vyr; wr; axr; ayr; alphar; ...]


    % Baumgarte gains
    if ~isfield(params, 'alpha_baum'), params.alpha_baum = 5.0;  end;
    if ~isfield(params, 'beta_baum'),  params.beta_baum  = 10.0; end;
    alpha = params.alpha_baum;
    beta  = params.beta_baum;

    % unpack state
    xr = x(1);  yr = x(2); thetar  = x(3);
    vxr = x(4); vyr = x(5); wr     = x(6);
    xt = x(7);  yt = x(8); thetat  = x(9);
    vxt = x(10); vyt = x(11); wt   = x(12);

    % get matrices
    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;        % 6x6
    B    = dyn.B;        % 6x2
    J    = dyn.J;        % 3x6
    dotJ = dyn.dotJ;     % 3x6
    vgen = dyn.v;        % 6x1

    % === position-level constraint h(q) ===
    % same geometry as in towing_dynamics_mats
    d  = params.d;
    Lt = params.Lt;

    r_rh_body = [-d;   0];
    r_th_body = [ Lt/2; 0];

    cr = cos(thetar); sr = sin(thetar);
    ct = cos(thetat); st = sin(thetat);

    Rr = [ cr, -sr;
           sr,  cr];
    Rt = [ ct, -st;
           st,  ct];

    p_r = Rr * r_rh_body;
    p_t = Rt * r_th_body;

    p_h_r  = [xr; yr] + p_r;
    p_h_tr = [xt; yt] + p_t;

    phi_pos = p_h_r - p_h_tr;    % 2x1 holonomic position error

    % third constraint (nonholonomic) has no position level: h3 = 0
    h = [phi_pos;
         0];                     % 3x1

    % === Baumgarte-stabilized constraint equation ===
    %   J a + dotJ v + 2 alpha J v + beta^2 h = 0
    Jv      = J * vgen;                           % 3x1
    rhs_con = -dotJ * vgen - 2*alpha*Jv - (beta^2) * h;  % 3x1

    % === KKT system ===
    % [ M   J';   [a    = [ B*u
    %   J   0 ]   lambda]   rhs_con ]
    m_con = size(J,1);
    K   = [M,  J.';
           J,  zeros(m_con, m_con)];

    rhs = [B*u;
           rhs_con];

    sol = K \ rhs;
    a   = sol(1:6);      % generalized acceleration

    % === assemble xdot ===
    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end
