clear; clc; close all;
addpath('D:\casADi') %% the path of your casadi
% Simulation parameters
tspan = 0:0.2:30;

% Geometry and dynamics params
params.d  = 0.134;
params.Lr = 0.50;  params.Wr = 0.19;
params.Lt = 0.514; params.Wt = 0.639;
params.m_r = 12;   params.I_r = 5;
params.m_t = 6.3;  params.I_t = 2;
params.Fmax= 1e3
params.Taumax= 1e3; 
% Initial state
params.x0 = zeros(12,1);   % [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
params.x0(1:3) = [0; 0; 0];  
xr0 = params.x0(1); 
yr0 = params.x0(2); 
thetar0 = params.x0(3);

% Robot hitch point in world coordinates
r_r_h = [xr0; yr0] - params.d*[cos(thetar0); sin(thetar0)];

% Trailer COM must satisfy:
%   r_th = r_r_h  = [xt; yt] + R(thetat)*[Lt/2 ; 0]
% Choose initial trailer heading same as robot
xt0 = r_r_h(1) - (params.Lt/2)*cos(thetar0);
yt0 = r_r_h(2) - (params.Lt/2)*sin(thetar0);

params.x0(7:9)   = [xt0; yt0; thetar0];  % initial trailer pose
params.x0(10:12) = [0;0;0];              % trailer initially at rest

% Use the feedback-linearizing controller via a wrapper that matches the simulator interface.
% controller = @(x,t) user_force_input(x,t);
%controller = @(x,t) fb_lin_wrapper(x, t, params);
%%initialize TO
dt_to  = tspan(2) - tspan(1);
N_to   = numel(tspan) - 1; 
y0=[-0.391;0]; %hitch's position
ref_to = make_straight_ref(tspan, y0);
save('ref_to.mat','ref_to');
%% Get Past TO result to reuse Warning: delete TO_Output.mat for new ref path
if ~isfile('TO_Output.mat')
    % load cached TO result
    S = load('TO_Output.mat', 'x_to', 'u_to', 'to_dbg');
    x_to   = S.x_to;
    u_to   = S.u_to;
    to_dbg = S.to_dbg;
else
    % run TO and save for next time
    [x_to, u_to, to_dbg] = towing_trajopt(dt_to, N_to, ref_to, params.x0, params);

    % ---- build time grid for to ----
    t_to = (0:N_to).' * dt_to;   % (N_to+1)x1
    
    % ---- normalize shapes to match build_to_log expectations ----
    X_to = x_to;
    if size(X_to,1) == 12 && size(X_to,2) == numel(t_to)
        X_to = X_to.';           % 12x(N+1) -> (N+1)x12
    end
    
    U_to = u_to;
    if ~isempty(U_to) && size(U_to,2) ~= 2 && size(U_to,1) == 2
        U_to = U_to.';           % 2xN -> Nx2
    end
    
    % ---- build to_log for analysis ----
    to_log = build_to_log(t_to, X_to, U_to, params);
    
    % ---- save everything ----
    save('TO_Output.mat', 'x_to', 'u_to', 'to_dbg', 't_to', 'to_log');
    end

%%
traj_to = build_to_qp_traj(dt_to, x_to, u_to, params);
params.to_traj = traj_to;
controller = @(x,t) qp_to_wrapper(x, t, params);

% Run simulation
simulate_planar_towing_full_dynamics(controller, tspan, params);
load('TO_Output.mat','to_log','to_dbg');
load('sim_log.mat','sim_log');
load('ref_to.mat','ref_to');
analyze_towing_results(sim_log, to_log, ref_to, params, to_dbg);
function to_log = build_to_log(t_to, X_to, U_to, params)
% build_to_log (minimal)
% packs TO outputs into the same schema used by analysis
%
% inputs:
%   t_to : (N x 1) or (1 x N)
%   X_to : (N x 12) or (12 x N)
%   U_to : (N-1 x 2) or (2 x (N-1)) or (N x 2) optional
%   params.d, params.Lt required

    t = t_to(:);
    N = numel(t);

    % normalize X to Nx12
    X = X_to;
    if size(X,1) == 12 && size(X,2) == N
        X = X.';  % 12xN -> Nx12
    end
    if size(X,1) ~= N || size(X,2) ~= 12
        error('X_to must be Nx12 or 12xN with N = numel(t_to)');
    end

    % normalize U to Nx2 (pad last sample if needed)
    U_pad = [];
    if nargin >= 3 && ~isempty(U_to)
        U = U_to;
        if size(U,2) ~= 2 && size(U,1) == 2
            U = U.'; % 2xK -> Kx2
        end

        if size(U,1) == N-1 && size(U,2) == 2
            U_pad = [U; U(end,:)];     % -> Nx2
        elseif size(U,1) == N && size(U,2) == 2
            U_pad = U;                 % already Nx2
        else
            error('U_to must be (N-1)x2, Nx2, or 2x(N-1)');
        end
    end

    d  = params.d;
    Lt = params.Lt;

    xr = X(:,1);  yr = X(:,2);  thetar = X(:,3);
    wr = X(:,6);
    xt = X(:,7);  yt = X(:,8);  thetat = X(:,9);

    p_r  = [xr, yr];
    p_tr = [xt, yt];

    p_h_r = [xr - d*cos(thetar),        yr - d*sin(thetar)];
    p_h_t = [xt + (Lt/2)*cos(thetat),   yt + (Lt/2)*sin(thetat)];

    e_h = sqrt(sum((p_h_r - p_h_t).^2, 2));

    to_log = struct();
    to_log.t     = t;
    to_log.X     = X;
    to_log.U     = U_pad;
    to_log.theta = thetar;
    to_log.wr    = wr;
    to_log.p_r   = p_r;
    to_log.p_tr  = p_tr;
    to_log.p_h_r = p_h_r;
    to_log.p_h_t = p_h_t;
    to_log.e_h   = e_h;
end
