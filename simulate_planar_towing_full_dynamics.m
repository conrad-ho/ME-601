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
dt_ctrl = 0.01;  % 100 Hz, might need to raise to 500Hz
[t, X, U, dbg] = simulate_sample_hold_ode45(controller_handle, tspan, params.x0, params, dt_ctrl, opts);

%%
% =========================================================
% build minimal analysis log (no for-loops)
% keeps: yaw, hitch error, paths, inputs, rv
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
sim_log.rv     = dbg.rv(:);               % [N x 1]  (speed-level constraint residual)

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
% - (optional) velocity projection each interval end (plan A)

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
    dbg.rv = zeros(numel(T),1);

    for k = 1:numel(T)-1
        tk = T(k);
        tk1 = T(k+1);

        % --- solve controller ONCE ---
        [F_drive, tau_r] = controller_handle(x, tk);
        u = [F_drive; tau_r];
        U(k,:) = u.';

        % --- integrate dynamics with u held constant ---
        f = @(tt,xx) towing_dynamic(xx, params, 'full').xdot(u);
        [~, xseg] = ode45(f, [tk tk1], x, opts);

        x = xseg(end,:).';

        % --- plan A: project velocity to satisfy Jv=0 ---
        x = project_velocity_to_constraints(x, params);

        X(k+1,:) = x.';

        % --- log rv ---
        dm = towing_dynamic(x, params, 'mats');
        dbg.rv(k+1) = norm(dm.J * dm.v);
    end
end
function x = project_velocity_to_constraints(x, params)
% project_velocity_to_constraints
% project generalized velocity v so that J*v = 0 (minimum correction in M-metric)
%
% state x: 12x1
%   [xr; yr; thetar; vxr; vyr; wr;  xt; yt; thetat; vxt; vyt; wt]
%
% assumes towing_dynamic(x, params, 'mats') returns:
%   dyn.M (6x6), dyn.J (3x6), dyn.v (6x1) = [vxr; vyr; wr; vxt; vyt; wt]

    dyn = towing_dynamic(x, params, 'mats');

    M = dyn.M;
    J = dyn.J;
    v = dyn.v;

    % solve MinvJt = M^{-1} J^T without forming inv(M)
    MinvJt = M \ (J.');          % 6x3

    % K = J M^{-1} J^T
    K = J * MinvJt;              % 3x3

    % velocity projection: v <- v - M^{-1}J^T K^{-1} (J v)
    v_corr = v - MinvJt * (K \ (J * v));

    % write back to x
    x(4:6)   = v_corr(1:3);
    x(10:12) = v_corr(4:6);
end
