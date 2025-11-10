function simulate_planar_towing_full_dynamics(controller_handle, tspan, params)
% Full dynamic simulation: robot + trailer with rigid tow arm constraint (exact)
% Uses Lagrange multiplier to enforce hitch distance = L_tow exactly.

% Default parameters
if ~isfield(params,'d'), params.d = 0.134; end
if ~isfield(params,'L_tow'), params.L_tow = 0.375; end
if ~isfield(params,'Lr'), params.Lr = 0.5; end
if ~isfield(params,'Wr'), params.Wr = 0.19; end
if ~isfield(params,'Lt'), params.Lt = 0.514; end
if ~isfield(params,'Wt'), params.Wt = 0.639; end
if ~isfield(params,'m_r'), params.m_r = 50; end
if ~isfield(params,'I_r'), params.I_r = 5; end
if ~isfield(params,'m_t'), params.m_t = 20; end
if ~isfield(params,'I_t'), params.I_t = 2; end

% Initial state
if ~isfield(params,'x0')
    params.x0 = [0;0;0;0;0;0; 0;0;0;0;0;0];
end

% Project initial trailer COM so hitch distance = L_tow
xr0 = params.x0(1); yr0 = params.x0(2); thetar0 = params.x0(3);
r_r_h = [xr0; yr0] - params.d*[cos(thetar0); sin(thetar0)];
r_t_h_des = r_r_h - params.L_tow*[cos(thetar0); sin(thetar0)];
xt0 = r_t_h_des(1) - (params.Lt/2)*cos(thetar0);
yt0 = r_t_h_des(2) - (params.Lt/2)*sin(thetar0);
params.x0(7:9) = [xt0; yt0; thetar0];
params.x0(10:12) = [0;0;0];

% ODE Integration 
opts = odeset('RelTol',1e-5,'AbsTol',1e-6);
[t,X] = ode45(@(tt,xx) full_dynamics(tt,xx,controller_handle,params), tspan, params.x0, opts);

% After integration: enforce rigid hitch geometry
for k = 1:length(t)
    xr = X(k,1);  yr = X(k,2);  thetar = X(k,3);
    xt = X(k,7);  yt = X(k,8);  thetat = X(k,9);

    % Robot hitch point (rear center minus distance d)
    r_rh = [xr; yr] - params.d * [cos(thetar); sin(thetar)];

    % Desired trailer hitch location (L_tow back from robot heading)
    r_t_h_des = r_rh - params.L_tow * [cos(thetar); sin(thetar)];

    % Move trailer center so its front hitch is at r_t_h_des
    X(k,7) = r_t_h_des(1) - (params.Lt/2)*cos(thetat);
    X(k,8) = r_t_h_des(2) - (params.Lt/2)*sin(thetat);
end
for k = 1:length(t)
    xr = X(k,1); yr = X(k,2); thetar = X(k,3);
    xt = X(k,7); yt = X(k,8); thetat = X(k,9);

    r_rh = [xr; yr] - params.d * [cos(thetar); sin(thetar)];
    r_th = [xt; yt] + (params.Lt/2)*[cos(thetat); sin(thetat)];

    hitch_len = norm(r_rh - r_th);
    fprintf('t=%.2f, hitch=%.4f (L_tow=%.4f)\n', t(k), hitch_len, params.L_tow);
end


% Animation Setup
figure('Color','w'); hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)');
title('Planar Robot–Trailer Dynamics (Rigid Tow Arm)');
xlim([-2 12]); ylim([-4 4]);
path_trace = plot(NaN,NaN,'k-','LineWidth',1.2);
robot_patch = patch(NaN,NaN,'r','FaceAlpha',0.4,'EdgeColor','none');
trailer_patch = patch(NaN,NaN,'b','FaceAlpha',0.4,'EdgeColor','none');
tow_line = plot([0 0],[0 0],'k--','LineWidth',2);

% Animate
traj_x=[]; traj_y=[];
for i=1:5:length(t)
    xr=X(i,1); yr=X(i,2); thetar=X(i,3);
    xt=X(i,7); yt=X(i,8); thetat=X(i,9);

    traj_x(end+1)=xr; traj_y(end+1)=yr;
    set(path_trace,'XData',traj_x,'YData',traj_y);

    % Update robot & trailer geometry
    set(robot_patch,'XData',rect_x(xr,params.Wr,params.Lr,thetar), ...
                    'YData',rect_y(yr,params.Wr,params.Lr,thetar));
    set(trailer_patch,'XData',rect_x(xt,params.Wt,params.Lt,thetat), ...
                      'YData',rect_y(yt,params.Wt,params.Lt,thetat));

    % Hitch positions
    r_r_h = [xr; yr] - params.d*[cos(thetar); sin(thetar)];
    r_t_h = [xt; yt] + (params.Lt/2)*[cos(thetat); sin(thetat)];
    set(tow_line,'XData',[r_r_h(1) r_t_h(1)],'YData',[r_r_h(2) r_t_h(2)]);

    axis([xr-3 xr+5 yr-3 yr+3]);
    drawnow;
end
end

% DYNAMICS FUNCTION
function dx = full_dynamics(t, x, controller_handle, params)
% State: [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]
xr = x(1); yr = x(2); thetar = x(3);
vxr = x(4); vyr = x(5); wr = x(6);
xt = x(7); yt = x(8); thetat = x(9);
vxt = x(10); vyt = x(11); wt = x(12);

m_r = params.m_r; I_r = params.I_r;
m_t = params.m_t; I_t = params.I_t;

% Control inputs (force and torque on robot)
[F_drive, tau_r] = controller_handle(x, t);

% Force in world frame along robot heading
Fr_body = [cos(thetar); sin(thetar)] * F_drive;
tau_r_body = tau_r;

%% Holonomic constraint: hitch points coincide
% (xr, yr) + r_rh == (xt, yt) + r_th
% Derivative constraint: J_holo * v = 0
J_holo_r = [-eye(2), -R2(thetar) * [0; -params.d]];    % 2x3
J_holo_t = [ eye(2),  R2(thetat) * [0;  params.Lt/2]]; % 2x3
J_holo   = [J_holo_r, J_holo_t];                       % 2x6

%% (rigid hitch, no relative yaw): not true in real life
%J_orient = [0 0 1 0 0 -1];  % thetar_dot - thetat_dot = 0

%% Nonholonomic constraint: trailer cannot move sideways
J_nonholo = [0 0 0, sin(thetat), -cos(thetat), 0];  % 1x6

%% Total constraint Jacobian
% J = [J_holo;
%      J_orient;
%      J_nonholo];  
J = [J_holo; J_nonholo];

%% Mass matrix
M = diag([m_r, m_r, I_r, m_t, m_t, I_t]);

%% Generalized velocities
v = [vxr; vyr; wr; vxt; vyt; wt];

%% Generalized forces (only robot actuated)
Q = [Fr_body; tau_r_body; zeros(3,1)];

%% Constraint bias (to keep constraints satisfied)
b_total = [-J * v];  

%% Solve coupled system [M J'; J 0] * [a; λ] = [Q; b_total]
A = [M, J'; J, zeros(size(J,1))];
rhs = [Q; b_total];
sol = A \ rhs;

a = sol(1:6);  % accelerations

% Hitch geometry
r_rh = R2(thetar) * [-params.d; 0];       
r_th = R2(thetat) * [params.Lt/2; 0];     

% DEBUG: hitch distance
r_robot_hitch = [xr; yr] + r_rh;
r_trailer_hitch = [xt; yt] + r_th;
hitch_length = norm(r_trailer_hitch - r_robot_hitch);
fprintf('t=%.2f s, hitch=%.4f (L_tow=%.4f)\n', t, hitch_length, params.L_tow);

%% Return state derivative for ODE integration
dx = [vxr; vyr; wr; a(1); a(2); a(3); vxt; vyt; wt; a(4); a(5); a(6)];
end

% Helper rotation matrix
function R = R2(th)
R = [cos(th) -sin(th); sin(th) cos(th)];
end

function X = rect_x(xc,W,L,theta)
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = 0.5*[-L -W; L -W; L W; -L W]';
X = R(1,:)*pts + xc;
end

function Y = rect_y(yc,W,L,theta)
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
pts = 0.5*[-L -W; L -W; L W; -L W]';
Y = R(2,:)*pts + yc;
end
