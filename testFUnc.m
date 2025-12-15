
tspan = 0:0.02:30;
clear functions   % 清所有 persistent（最狠但最有效）
% Geometry and dynamics params
params.d  = 0.134;
params.Lr = 0.50;  params.Wr = 0.19;
params.Lt = 0.514; params.Wt = 0.639;
params.m_r = 12;   params.I_r = 5;
params.m_t = 6.3;  params.I_t = 2;

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

params.to_traj = make_hold_to_traj(tspan, params.x0, params);
controller = @(x,t) qp_to_wrapper(x, t, params);

% Run simulation
simulate_planar_towing_full_dynamics(controller, tspan, params);

function traj = make_hold_to_traj(tspan, x0, params)
% make a placeholder TO trajectory that holds the initial hitch point
% format matches params.to_traj:
%   t: 1x(N+1)
%   x: (N+1)x12
%   u: Nx2
%   y, ydot, yddot_ff: 2x(N+1)

    t = tspan(:).';          % 1x(N+1)
    N = numel(t) - 1;

    % constant state guess (not actually used if wrapper only reads y/ydot/yddot_ff)
    x = repmat(x0(:).', N+1, 1);    % (N+1)x12
    u = zeros(N, 2);                % Nx2

    % compute hitch at initial state
    dyn0 = towing_dynamic(x0, params, 'mats');
    y0 = dyn0.y;                    % 2x1

    y = repmat(y0, 1, N+1);         % 2x(N+1)
    ydot = zeros(2, N+1);
    yddot_ff = zeros(2, N+1);

    traj = struct;
    traj.t = t;
    traj.x = x;
    traj.u = u;
    traj.y = y;
    traj.ydot = ydot;
    traj.yddot_ff = yddot_ff;
end
