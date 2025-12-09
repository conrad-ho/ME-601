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
params.alpha_baum = 5.0;
params.beta_baum  = 10.0;
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
p_tr0=[xt0;yt0];

opts.theta0 = 0;          % rad
opts.d  = 0.134;          % robot COM -> hitch
opts.Lt = 0.257;          % trailer COM -> hitch (= Lt/2)

ref_guess = make_to_ref_Sarc_robot_trailer(tspan, p_tr0, opts);
save('log.mat','ref_guess'); 
%% Get Past TO result to reuse Warning: delete TO_Output.mat for new ref path

    % run TO and save for next time
    x0_full = params.x0;   % 12x1, 原来的 full state
    x0_kin  = [x0_full(1);   % xr
           x0_full(2);   % yr
           x0_full(3);   % thetar
           x0_full(9)];  % theta_t (trailer heading)
    [x_kin, u_kin, to_dbg] = towing_trajopt_kinematic(dt_to, N_to, ref_guess, x0_kin, params);
    ref_to = lift_kin_to_full_ref(x_kin, params, dt_to);% 3) 用 TO 结果 lift 出最终要跟踪的 ref_to
    %check_to_rigid_constraint(x_kin, params);   % 这是测 4 维版本
    x_full_ref = lift_kin_to_full_state(x_kin, params);
    %check_to_rigid_constraint(x_full_ref, params);
   %visualize_to_vs_ref(x_full_ref, ref_guess, params, dt_to);
    %% QP
    traj_to        = build_to_qp_traj(dt_to, x_full_ref, u_kin, params);
    params.to_traj = traj_to;
    params.ref_to=ref_to;
    params.u_min=[-30,-5];
    params.u_max=[30,5];
    params.u_const = [0; 0];   
    %controller     = @(x,t) qp_to_wrapper(x, t, params);
    params.qp_verbose = true;      % 或调成 false 静音
params.dt_qp      = 0.02;      % 50 Hz 更新一次控制

% 仿真设置
params.dt_sim = 1e-3;         % 比如 1ms
params.u_const = [0; 0];      % 零输入

controller = @(x,t) qp_path_tracker_from_to_safe(x, t, params);

log = simulate_planar_towing_full_dynamics_fixed(controller, [0, 30], params);



% 后处理，例如：
t   = log.t;
X   = log.X;
u   = log.u;
    
 %%
function [F_drive, tau_r, dbg] = const_controller(~, ~, params)
%CONST_CONTROLLER
%   Very simple controller: constant inputs (from params or default)

    if isfield(params,'u_const')
        u = params.u_const(:);
        if numel(u) ~= 2
            error('const_controller:bad_u_const', ...
                  'params.u_const must be 2x1.');
        end
        F_drive = u(1);
        tau_r   = u(2);
    else
        F_drive = 0;
        tau_r   = 0;
    end

    dbg = struct();   % dummy
end

%%
   function x_full = lift_kin_to_full_state(x_kin, params)
%LIFT_KIN_TO_FULL_STATE
%   Lift kinematic trajectory [xr, yr, thetar, theta_t]
%   into a fake full-state trajectory (N+1) x 12 for diagnostics:
%
%   x_full(k,:) = [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
%
%   where:
%     - trailer COM (xt,yt,thetat) is reconstructed from hitch geometry:
%         p_h  = [xr - d*cos(thetar);
%                 yr - d*sin(thetar)];
%         p_tr = p_h - [Lt*cos(theta_t);
%                       Lt*sin(theta_t)];
%
%   velocities are set to zeros (this is OK for checking rigid constraint).

    if size(x_kin,2) ~= 4
        error('x_kin must be (N+1) x 4 = [xr, yr, thetar, theta_t].');
    end
    if ~isfield(params,'d') || ~isfield(params,'Lt')
        error('params must contain d and Lt.');
    end

    d  = params.d;
    Lt = params.Lt;

    Np = size(x_kin,1);  % N+1

    x_full = zeros(Np, 12);

    for k = 1:Np
        xr      = x_kin(k,1);
        yr      = x_kin(k,2);
        thetar  = x_kin(k,3);
        theta_t = x_kin(k,4);

        % robot COM
        xtmp = xr;
        ytmp = yr;

        % hitch from robot
        p_h = [xr - d*cos(thetar);
               yr - d*sin(thetar)];

        % trailer COM (left of hitch)
        p_tr = p_h - [Lt*cos(theta_t);
                      Lt*sin(theta_t)];

        xt = p_tr(1);
        yt = p_tr(2);
        thetat = theta_t;

        % pack into 12-dim state: [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
        x_full(k,:) = [xr, yr, thetar, ...
                       0,  0,  0, ...      % vxr, vyr, wr (unused for rigid check)
                       xt, yt, thetat, ...
                       0,  0,  0];         % vxt, vyt, wt
    end
end
