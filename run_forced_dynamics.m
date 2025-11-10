clear; clc; close all;

% Simulation parameters
tspan = 0:0.02:30;

% Geometry and dynamics params
params.d = 0.134; params.L_tow = 0.375;
params.Lr = 0.50; params.Wr = 0.19;
params.Lt = 0.514; params.Wt = 0.639;
params.m_r = 12; params.I_r = 5;
params.m_t = 6.3; params.I_t = 3;

% Initial state
params.x0 = zeros(12,1); % [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
params.x0(1:3) = [0; 0; 0]; 

xr0 = params.x0(1); yr0 = params.x0(2); thetar0 = params.x0(3);
r_r_h = [xr0; yr0] - params.d*[cos(thetar0); sin(thetar0)];
r_t_h_des = r_r_h - params.L_tow*[cos(thetar0); sin(thetar0)];
xt0 = r_t_h_des(1) - (params.Lt/2)*cos(thetar0);
yt0 = r_t_h_des(2) - (params.Lt/2)*sin(thetar0);
params.x0(7:9) = [xt0; yt0; thetar0]; % trailer initial pose
params.x0(10:12) = [0;0;0]; % trailer initially still
controller = @(x,t) user_force_input(x,t);
simulate_planar_towing_full_dynamics(controller, tspan, params);

%% Force/Torque Input Function
function [F, tau] = user_force_input(~,t)
    % test: constant thrust, small steering pulse
    F = 30;   % [N] forward drive
    if t > 4 && t < 6
        tau = 1;  % small torque to turn
    else
        tau = 0;
    end
end
