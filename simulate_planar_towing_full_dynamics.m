function simulate_planar_towing_full_dynamics(controller_handle, tspan, params)
% full dynamic simulation: robot + trailer with rigid hitch constraint.
% dynamics is provided by towing_dynamic(x, params, 'full').

if ~isfield(params,'d'),  params.d  = 0.134; end
if ~isfield(params,'Lr'), params.Lr = 0.5;   end
if ~isfield(params,'Wr'), params.Wr = 0.19;  end
if ~isfield(params,'Lt'), params.Lt = 0.514; end
if ~isfield(params,'Wt'), params.Wt = 0.639; end
if ~isfield(params,'m_r'), params.m_r = 50;  end
if ~isfield(params,'I_r'), params.I_r = 5;   end
if ~isfield(params,'m_t'), params.m_t = 20;  end
if ~isfield(params,'I_t'), params.I_t = 2;   end

% initial state default
if ~isfield(params,'x0')
    params.x0 = zeros(12,1);
end

% project trailer initial state so the hitch points coincide exactly
xr0 = params.x0(1); yr0 = params.x0(2); thetar0 = params.x0(3);

% robot hitch in world
r_r_h = [xr0; yr0] - params.d*[cos(thetar0); sin(thetar0)];

% trailer com must satisfy: r_th = r_r_h, with thetat0 = thetar0
xt0 = r_r_h(1) - (params.Lt/2)*cos(thetar0);
yt0 = r_r_h(2) - (params.Lt/2)*sin(thetar0);
params.x0(7:9)   = [xt0; yt0; thetar0];
params.x0(10:12) = [0;0;0];

% integrate
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
dt_ctrl = 0.01;  % 100 Hz

[t, X, U, dbg] = simulate_sample_hold_ode45(controller_handle, tspan, params.x0, params, dt_ctrl, opts);

%%
% =========================================================
% build minimal analysis log (no for-loops)
% keeps: yaw, hitch error, paths, inputs, rv/qp_rv/qp_ra
% =========================================================

N  = numel(t);

% state slices
xr     = X(:,1);   yr     = X(:,2);   thetar = X(:,3);
wr     = X(:,6);
xt     = X(:,7);   yt     = X(:,8);   thetat = X(:,9);

d  = params.d;
Lt = params.Lt;

% paths
p_r  = [xr, yr];                      % robot com
p_tr = [xt, yt];                      % trailer com

% hitch points (robot hitch and trailer hitch)
p_h_r = [xr - d*cos(thetar),  yr - d*sin(thetar)];                 % robot-side hitch point
p_h_t = [xt + (Lt/2)*cos(thetat), yt + (Lt/2)*sin(thetat)];        % trailer-side hitch point

% hitch position error (should be ~0)
e_h = sqrt(sum((p_h_r - p_h_t).^2, 2));    % [N x 1]

% inputs: U is (N-1)x2 from sample/hold; pad to N for plotting convenience
if isempty(U)
    U_pad = zeros(N,2);
else
    U_pad = [U; U(end,:)];                % [N x 2]
end

% pack log
sim_log = struct();
sim_log.t      = t(:);                    % [N x 1]
sim_log.X      = X;                       % [N x 12]
sim_log.U      = U_pad;                   % [N x 2]  (F_drive, tau_r)
sim_log.theta  = thetar;                  % [N x 1]
sim_log.wr     = wr;                      % [N x 1]
sim_log.p_r    = p_r;                     % [N x 2]
sim_log.p_tr   = p_tr;                    % [N x 2]
sim_log.p_h_r  = p_h_r;                   % [N x 2]
sim_log.p_h_t  = p_h_t;                   % [N x 2]
sim_log.e_h    = e_h;                     % [N x 1]


% also store bounds for analysis convenience
if isfield(params,'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
if isfield(params,'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end
sim_log.lb = [-Fmax; -Taumax];
sim_log.ub = [ Fmax;  Taumax];
% constraint monitors
if isfield(dbg,'rv'), sim_log.rv = dbg.rv(:); else, sim_log.rv = nan(N,1); end

% qp monitors (numeric-only, safe to save)
if isfield(dbg,'qp_rv'), sim_log.qp_rv = dbg.qp_rv(:); else, sim_log.qp_rv = nan(N,1); end
if isfield(dbg,'qp_ra'), sim_log.qp_ra = dbg.qp_ra(:); else, sim_log.qp_ra = nan(N,1); end
if isfield(dbg,'qp_exitflag'), sim_log.qp_exitflag = dbg.qp_exitflag(:); else, sim_log.qp_exitflag = nan(N,1); end
if isfield(dbg,'qp_condK'), sim_log.qp_condK = dbg.qp_condK(:); else, sim_log.qp_condK = nan(N,1); end
% qp hitch monitors (vector signals)
if isfield(dbg,'qp_h'),       sim_log.qp_h       = dbg.qp_h;        end
if isfield(dbg,'qp_hdot'),    sim_log.qp_hdot    = dbg.qp_hdot;     end
if isfield(dbg,'qp_e_yddot'), sim_log.qp_e_yddot = dbg.qp_e_yddot;  end
if isfield(dbg,'qp_yddot'),   sim_log.qp_yddot   = dbg.qp_yddot;    end
if isfield(dbg,'u'),        sim_log.u        = dbg.u;        end
if isfield(dbg,'u_sat'),    sim_log.u_sat    = dbg.u_sat;    end
if isfield(dbg,'u_margin'), sim_log.u_margin = dbg.u_margin; end
if isfield(dbg,'qp_Ay'), sim_log.A_y = dbg.qp_Ay; end

save('sim_log.mat','sim_log');

%%
animate_planar_towing_traj(t, X, params, opts)
end

% =========================================================
% ODE RHS using unified towing_dynamic interface
% =========================================================
function [T, X, U, dbg] = simulate_sample_hold_ode45(controller_handle, tspan, x0, params, dt_ctrl, opts)
% simulate with sample-and-hold control:
% - solve controller once every dt_ctrl
% - integrate dynamics with u held constant within each interval using ode45
% - project velocity at interval end so that Jv=0

if nargin < 6 || isempty(opts)
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
end

t0 = tspan(1);
tf = tspan(end);

% build control grid
N = ceil((tf - t0)/dt_ctrl);
t_grid = t0 + (0:N)*dt_ctrl;
t_grid(end) = tf;   % ensure ends exactly at tf

x = x0(:);

T = t_grid(:);
X = zeros(numel(T), numel(x));
U = zeros(numel(T)-1, 2);

X(1,:) = x.';

dbg = struct;
dbg.rv          = zeros(numel(T),1);

% qp debug (store only doubles; dummy controller -> stays NaN)
dbg.qp_rv       = nan(numel(T),1);
dbg.qp_ra       = nan(numel(T),1);
dbg.qp_exitflag = nan(numel(T),1);
dbg.qp_condK    = nan(numel(T),1);

dbg.qp_h        = nan(numel(T),2);   % hitch pos error vector (y-yd)
dbg.qp_hdot     = nan(numel(T),2);   % hitch vel error vector (ydot-ydot_d)
dbg.qp_e_yddot  = nan(numel(T),2);   % accel error vector (yddot - yddot_des)  (preferred)

dbg.u        = nan(numel(T),2);    % u actually applied at each grid time (pad)
dbg.u_sat    = false(numel(T),2);  % [sat_F, sat_tau]
dbg.u_margin = nan(numel(T),2);    % min distance to bounds for each inp
dbg.qp_Ay = nan(numel(T), 2, 2);   % store A_y(t)

for k = 1:numel(T)-1
    tk  = T(k);
    tk1 = T(k+1);

    % --- solve controller ONCE (support optional qp_dbg) ---
    F_drive = 0; tau_r = 0; qp_dbg = [];  % safe defaults every step
    %try
        % try 3-output signature first
        [F_drive, tau_r, qp_dbg] = controller_handle(x, tk);
        %{
    catch
        try
            % fallback: 2-output signature
            [F_drive, tau_r] = controller_handle(x, tk);
            qp_dbg = [];
        catch
            % keep defaults if controller fails
            qp_dbg = [];
        end
    end
        %}
    u = [F_drive; tau_r];
    U(k,:) = u.';
    if ~isempty(qp_dbg) && isstruct(qp_dbg) && isfield(qp_dbg,'lb') && isfield(qp_dbg,'ub') ...
        && ~isempty(qp_dbg.lb) && ~isempty(qp_dbg.ub)
    lb = double(qp_dbg.lb(:));
    ub = double(qp_dbg.ub(:));
    else
    if isfield(params,'Fmax'),   Fmax   = params.Fmax;   else, Fmax   = 1e3; end
    if isfield(params,'Taumax'), Taumax = params.Taumax; else, Taumax = 1e3; end
    lb = [-Fmax; -Taumax];
    ub = [ Fmax;  Taumax];
    end
    dbg.u(k+1,:) = u.';  % keep your k+1 alignment convention

    % use relative threshold (more robust than exact equality)
    Fmax_eff   = ub(1);
    Taumax_eff = ub(2);
    
    dbg.u_sat(k+1,1) = abs(u(1)) >= 0.995*Fmax_eff;
    dbg.u_sat(k+1,2) = abs(u(2)) >= 0.995*Taumax_eff;
    
    % distance to nearest bound (0 means exactly on bound)
    dbg.u_margin(k+1,:) = min((u - lb).', (ub - u).');
    % --- extract qp rv/ra (numeric only; NaN if unavailable) ---
    [qp_rv_k, qp_ra_k, qp_exitflag_k, qp_condK_k] = extract_qp_residuals(x, u, params, qp_dbg);
    dbg.qp_rv(k+1)       = qp_rv_k;
    dbg.qp_ra(k+1)       = qp_ra_k;
    dbg.qp_exitflag(k+1) = qp_exitflag_k;
    dbg.qp_condK(k+1)    = qp_condK_k;
    
    if ~isempty(qp_dbg) && isstruct(qp_dbg)
    if isfield(qp_dbg,'h') && ~isempty(qp_dbg.h)
        try, dbg.qp_h(k+1,:) = double(full(qp_dbg.h(:))).'; end
    end
    if isfield(qp_dbg,'hdot') && ~isempty(qp_dbg.hdot)
        try, dbg.qp_hdot(k+1,:) = double(full(qp_dbg.hdot(:))).'; end
    end

    % prefer accel error if you stored it as qp_dbg.e_yddot in controller
    if isfield(qp_dbg,'e_yddot') && ~isempty(qp_dbg.e_yddot)
        try, dbg.qp_e_yddot(k+1,:) = double(full(qp_dbg.e_yddot(:))).'; end
    end

    % optional: store actual yddot
    if isfield(qp_dbg,'yddot_opt') && ~isempty(qp_dbg.yddot_opt)
        try, dbg.qp_yddot(k+1,:) = double(full(qp_dbg.yddot_opt(:))).'; end
    end
    end
    if ~isempty(qp_dbg) && isstruct(qp_dbg)
    if isfield(qp_dbg,'A_y') && ~isempty(qp_dbg.A_y)
        dbg.qp_Ay(k+1,:,:) = double(full(qp_dbg.A_y));
    end
    end
    % --- integrate dynamics with u held constant ---
    dyn_full = towing_dynamic(x, params, 'full');
    f = @(tt,xx) dyn_full.xdot(u);
    [~, xseg] = ode45(f, [tk tk1], x, opts);
    x = xseg(end,:).';
    
    % --- project position/orientation to satisfy g(q)=0 (rigid hitch) ---
    x = project_position_to_constraints(x, params);

    % --- project velocity to satisfy Jv=0 ---
    x = project_velocity_to_constraints(x, params);
    % hard fail if NaN (do not hide bugs)
    if any(~isfinite(x))
    error('NaN/Inf after projection at step k=%d, t=%.6f', k, tk1);
    end
    X(k+1,:) = x.';

    % --- log rv from current state ---
    dm = towing_dynamic(x, params, 'mats');
    dbg.rv(k+1) = double(norm(full(dm.J * dm.v)));

    %=== testing, remove at will
    d  = params.d;
    xr = x(1);
    yr = x(2);
    th = x(3);

    p_hr = [ xr - d*cos(th);
             yr - d*sin(th) ];
    Lt = params.Lt;
    xt = x(7);
    yt = x(8);
    th = x(9);
    %{
    p_ht = [ xt + (Lt/2)*cos(th);
             yt + (Lt/2)*sin(th) ];
    e_h  = norm(p_hr - p_ht);
    save('log.mat','e_h')
    %}
end
end

function x = project_velocity_to_constraints(x, params)
% project generalized velocity v so that J*v = 0 (minimum correction in M-metric)

dyn = towing_dynamic(x, params, 'mats');

M = dyn.M;
J = dyn.J;
v = dyn.v;

MinvJt = M \ (J.');          % 6x3
K = J * MinvJt;              % 3x3

v_corr = v - MinvJt * (K \ (J * v));

x(4:6)   = v_corr(1:3);
x(10:12) = v_corr(4:6);
end

function x = project_position_to_constraints(x, params)
% project positions onto rigid hitch constraint manifold
% enforces:
%   1) hitch points coincide in x/y
%   2) trailer yaw = robot yaw (rigid hitch orientation)

% unpack
d   = params.d;
Lt  = params.Lt;

xr  = x(1);  yr  = x(2);  thetar = x(3);

% robot hitch position in world
p_hr = [xr - d*cos(thetar);
        yr - d*sin(thetar)];

% enforce rigid hitch orientation: thetat = thetar
x(9) = thetar;

% set trailer com so that trailer hitch matches robot hitch exactly
% trailer hitch: p_ht = [xt; yt] + (Lt/2)*[cos(thetat); sin(thetat)]
xt = p_hr(1) - (Lt/2)*cos(thetar);
yt = p_hr(2) - (Lt/2)*sin(thetar);

x(7) = xt;
x(8) = yt;

% fail fast on bad numbers
if any(~isfinite(x([1:3,7:9])))
    error('non-finite state in position projection');
end
end

function [rv, ra, exitflag, condK] = extract_qp_residuals(x, u, params, qp_dbg)
% return scalar doubles; NaN if not available
rv = NaN; ra = NaN; exitflag = NaN; condK = NaN;

% 0) if qp_dbg exists, try to read directly
if ~isempty(qp_dbg) && isstruct(qp_dbg)

    if isfield(qp_dbg,'exitflag') && ~isempty(qp_dbg.exitflag)
        exitflag = double(qp_dbg.exitflag);
    end
    if isfield(qp_dbg,'condK') && ~isempty(qp_dbg.condK)
        condK = double(qp_dbg.condK);
    end

    % prefer residual vectors if provided
    if isfield(qp_dbg,'res_vel') && ~isempty(qp_dbg.res_vel)
        try
            rv = double(norm(full(qp_dbg.res_vel)));
        catch
        end
    end
    if isfield(qp_dbg,'res_acc') && ~isempty(qp_dbg.res_acc)
        try
            ra = double(norm(full(qp_dbg.res_acc)));
        catch
        end
    end

    % sometimes store scalars directly
    if isnan(rv) && isfield(qp_dbg,'rv') && ~isempty(qp_dbg.rv)
        try
            rv = double(full(qp_dbg.rv));
        catch
        end
    end
    if isnan(ra) && isfield(qp_dbg,'ra') && ~isempty(qp_dbg.ra)
        try
            ra = double(full(qp_dbg.ra));
        catch
        end
    end
end

% 1) if rv still NaN, compute from mats (depends only on state)
if isnan(rv)
    dm = towing_dynamic(x, params, 'mats');
    rv = double(norm(full(dm.J * dm.v)));
end

% 2) if ra still NaN, compute accel-level residual using projected accel a(u)
% ra := || J*a + dotJ*v ||
if isnan(ra)
    try
        dp = towing_dynamic(x, params, 'proj');  % must provide J, dotJ, v, and a(u)
        a  = dp.a(u);
        ra = double(norm(full(dp.J * a + dp.dotJ * dp.v)));
    catch
        ra = NaN;
    end
end



end
