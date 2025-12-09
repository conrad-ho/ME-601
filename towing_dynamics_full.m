function xdot = towing_dynamics_full(x, u, params)
%TOWING_DYNAMICS_FULL
%   Continuous-time rigid robot+trailer dynamics with Baumgarte-stabilized
%   holonomic (hitch) + nonholonomic constraints.
%
% State:
%   x = [xr; yr; thetar; vxr; vyr; wr;
%        xt; yt; thetat; vxt; vyt; wt]  (12x1)
%
% Input:
%   u = [F_drive; tau_r] (2x1)
%
% Output:
%   xdot, same ordering as x:
%     [vxr; vyr; wr; axr; ayr; alphar; ...] (12x1)
%
% Debug options:
%   - if params.debug_full == true, extra info will be saved when NaN / near-singular KKT appears.

    % ---------- basic input checks ----------
    x = x(:);
    u = u(:);

    if numel(x) ~= 12
        error('towing_dynamics_full: state x must be 12x1, got %dx1', numel(x));
    end
    if numel(u) ~= 2
        error('towing_dynamics_full: input u must be 2x1, got %dx1', numel(u));
    end

    if any(~isfinite(x))
        error('towing_dynamics_full: input state x contains NaN/Inf');
    end
    if any(~isfinite(u))
        error('towing_dynamics_full: input control u contains NaN/Inf');
    end

    debug_on = isfield(params,'debug_full') && params.debug_full;

    % ---------- Baumgarte gains ----------
    if ~isfield(params, 'alpha_baum'), params.alpha_baum = 5.0;  end
    if ~isfield(params, 'beta_baum'),  params.beta_baum  = 10.0; end
    alpha = params.alpha_baum;
    beta  = params.beta_baum;

    % ---------- geometry parameters ----------
    if ~isfield(params,'d'),  error('params.d (hitch offset) missing');  end
    if ~isfield(params,'Lt'), error('params.Lt (trailer length) missing'); end

    d  = params.d;
    Lt = params.Lt;

    % ---------- unpack state ----------
    xr   = x(1);  yr   = x(2);  thetar = x(3);
    vxr  = x(4);  vyr  = x(5);  wr     = x(6);
    xt   = x(7);  yt   = x(8);  thetat = x(9);
    vxt  = x(10); vyt  = x(11); wt     = x(12); %#ok<NASGU>

    % ---------- dynamic matrices ----------
    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;        % 6x6
    B    = dyn.B;        % 6x2
    J    = dyn.J;        % 3x6
    dotJ = dyn.dotJ;     % 3x6
    vgen = dyn.v;        % 6x1 generalized velocity
    Hgen = dyn.H;        % 6x1 Coriolis/grav/etc

    % sanity check on dyn
    if any(~isfinite(M(:))) || any(~isfinite(B(:))) || ...
       any(~isfinite(J(:))) || any(~isfinite(dotJ(:))) || ...
       any(~isfinite(vgen(:))) || any(~isfinite(Hgen(:)))
        if debug_on
            dbg.x    = x;
            dbg.u    = u;
            dbg.M    = M;
            dbg.B    = B;
            dbg.J    = J;
            dbg.dotJ = dotJ;
            dbg.vgen = vgen;
            dbg.Hgen = Hgen;
            save('debug_dynamics_mats_nan.mat','dbg');
        end
        error('towing_dynamics_full: towing_dynamics_mats returned NaN/Inf');
    end

    % ---------- position-level hitch constraint h(q) ----------
    % same geometry as in towing_dynamics_mats
    r_rh_body = [-d;    0];      % robot body -> hitch (robot side)
    r_th_body = [ Lt/2; 0];      % trailer body -> hitch (trailer side)

    cr = cos(thetar); sr = sin(thetar);
    ct = cos(thetat); st = sin(thetat);

    Rr = [ cr, -sr;
           sr,  cr];
    Rt = [ ct, -st;
           st,  ct];

    p_r  = Rr * r_rh_body;       % hitch position in world, from robot
    p_t  = Rt * r_th_body;       % hitch position in world, from trailer

    p_h_r  = [xr; yr] + p_r;
    p_h_tr = [xt; yt] + p_t;

    phi_pos = p_h_r - p_h_tr;    % 2x1 holonomic position error

    % third constraint (nonholonomic) has no position level
    h = [phi_pos;
         0];                     % 3x1

    % ---------- Baumgarte-stabilized constraint equation ----------
    %   J a + dotJ v + 2 alpha J v + beta^2 h = 0
    Jv      = J * vgen;                           % 3x1
    rhs_con = -dotJ * vgen - 2*alpha*Jv - (beta^2) * h;  % 3x1

    if any(~isfinite(rhs_con))
        if debug_on
            dbg.x      = x;
            dbg.u      = u;
            dbg.J      = J;
            dbg.dotJ   = dotJ;
            dbg.vgen   = vgen;
            dbg.h      = h;
            dbg.rhs_con= rhs_con;
            save('debug_rhs_con_nan.mat','dbg');
        end
        error('towing_dynamics_full: rhs_con contains NaN/Inf');
    end

    % ---------- KKT system ----------
    % [ M   J';   [a    = [ B*u
    %   J   0 ]   lambda]   rhs_con ]
    m_con = size(J,1);
    K   = [M,  J.';
           J,  zeros(m_con, m_con)];

    rhs = [B*u;
           rhs_con];

    if any(~isfinite(K(:))) || any(~isfinite(rhs(:)))
        if debug_on
            dbg.x    = x;
            dbg.u    = u;
            dbg.M    = M;
            dbg.J    = J;
            dbg.K    = K;
            dbg.rhs  = rhs;
            save('debug_K_nan.mat','dbg');
        end
        error('towing_dynamics_full: K or rhs contains NaN/Inf before solve');
    end

    % ---------- check conditioning of K, regularize if needed ----------
    rc = rcond(K);

    if isnan(rc) || rc < 1e-10
        if debug_on
            warning('towing_dynamics_full: KKT matrix near singular, rcond = %g (saving debug_K_singular.mat)', rc);
            dbg.x    = x;
            dbg.u    = u;
            dbg.M    = M;
            dbg.J    = J;
            dbg.dotJ = dotJ;
            dbg.vgen = vgen;
            dbg.h    = h;
            dbg.rhs_con = rhs_con;
            dbg.K    = K;
            dbg.rhs  = rhs;
            dbg.rcond = rc;
            save('debug_K_singular.mat','dbg');
        else
            warning('towing_dynamics_full: KKT matrix near singular, rcond = %g', rc);
        end

        % small Tikhonov regularization to keep simulation from blowing up
        epsK = 1e-6;
        Kreg = K + epsK * eye(size(K));
        sol  = Kreg \ rhs;
    else
        sol = K \ rhs;
    end

    if any(~isfinite(sol))
        if debug_on
            dbg.x    = x;
            dbg.u    = u;
            dbg.K    = K;
            dbg.rhs  = rhs;
            dbg.sol  = sol;
            save('debug_sol_nan.mat','dbg');
        end
        error('towing_dynamics_full: solution of KKT system contains NaN/Inf');
    end

    a   = sol(1:6);      % generalized acceleration

    % ---------- assemble xdot ----------
    qdot = vgen;
    vdot = a;

    xdot = [qdot;
            vdot];
end
