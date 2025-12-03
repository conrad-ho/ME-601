function varargout = fb_lin_wrapper(x, t, params)
% Simple wrapper to call fb_lin_controller with default/desired trajectories.
% Returns [F_drive, tau_r] to match simulate_planar_towing_full_dynamics

% Trajectory 1: Hold current hitch location (stabilization)
% xr = x(1); yr = x(2); thetar = x(3);
% d = params.d;
% y0 = [ xr - d*cos(thetar); yr - d*sin(thetar) ];
% yd = y0;
% yd_dot = [0;0];
% yd_ddot = [0;0];

% Trajectory 2: follow a circular trajectory in world coords
% R = 1.0; omega = 0.05;
% center = [1.5; 0];
% yd = center + R*[cos(omega*t); sin(omega*t)];
% yd_dot = R*omega*[-sin(omega*t); cos(omega*t)];
% yd_ddot = -R*omega^2*[cos(omega*t); sin(omega*t)];

R = 1;
omega = 0.1;

% Get initial hitch point
xr = x(1); yr = x(2); thetar = x(3);
d = params.d;
y0 = [xr - d*cos(thetar); yr - d*sin(thetar)];

% Generate circle around current hitch position
yd = y0 + R*[cos(omega*t); sin(omega*t)];
yd_dot = R*omega*[-sin(omega*t); cos(omega*t)];
yd_ddot = -R*omega^2*[cos(omega*t); sin(omega*t)];


% Trajectory 3: straight-line desired hitch trajectory
% vx = 1;
% yd = [vx*t; 0];
% yd_dot = [vx; 0];
% yd_ddot = [0; 0];

% Call the feedback-linearizing controller
[F_drive, tau_r] = fb_lin_controller(x, t, params, yd, yd_dot, yd_ddot);

% Return outputs
varargout{1} = F_drive;
varargout{2} = tau_r;
end
