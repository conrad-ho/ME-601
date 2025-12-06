% experiment_spin_near_goal.m
clear; clc; close all;

%% 加载上一次完整仿真的 log
data = load('log.mat');   % 里面应有 log 结构
log_full = data.log;

% 终点 index（你也可以选倒数几步之一）
idx_end = size(log_full.X,1);

% 终点全状态 x_term (12x1)
x_term = log_full.X(idx_end,:).';

% 对应的 hitch 终点位置 (2x1)
y_goal = log_full.y(idx_end,:).';

%% 在终点状态上加一个 yaw 偏角
yaw_offset_deg = 45;                 % 你可以改成 20, 60 等等
yaw_offset_rad = deg2rad(yaw_offset_deg);

x0_pert = x_term;
x0_pert(3) = x0_pert(3) + yaw_offset_rad;   % 只改 robot 车头角 thetar

%% 设置仿真参数
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

params.x0 = x0_pert;        % 新仿真的初始状态

dt      = 0.02;
T_short = 6.0;              % 在终点附近再跑 6 秒
tspan   = 0:dt:T_short;

%% 定义一个 “把 hitch 固定在 y_goal 附近” 的 controller wrapper
% 这里假设你有 hitch_ref 和 task_space_qp_controller 两个函数
controller_handle = @(x,t) hold_hitch_at_goal_controller(x,t,params,y_goal);

%% 跑仿真
simulate_planar_towing_full_dynamics(controller_handle, tspan, params);

% 你可以在 simulate 里面再做日志/画图，或在这里额外 plot 一次 theta_r(t) 看有没有绕圈
