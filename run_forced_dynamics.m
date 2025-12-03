clear; clc; close all;

% Simulation parameters
tspan = 0:0.02:30;

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

% Use the feedback-linearizing controller via a wrapper that matches the simulator interface.
% controller = @(x,t) user_force_input(x,t);
controller = @(x,t) fb_lin_wrapper(x, t, params);


% Run simulation
simulate_planar_towing_full_dynamics(controller, tspan, params);

% Force/Torque Input Function, not being used
function [F, tau] = user_force_input(~,t)
    % test: constant thrust, small steering pulse
    F = 30;   % constant forward force
    if t > 4 && t < 6
        tau = 1;  % steering torque pulse
    else
        tau = 0;
    end
end