function simulate_planar_towing_full_dynamics(controller_handle, tspan, params)
% Full dynamic simulation: robot + trailer with rigid hitch constraint (exact).

if ~isfield(params,'d'),  params.d  = 0.134; end
if ~isfield(params,'Lr'), params.Lr = 0.5;   end
if ~isfield(params,'Wr'), params.Wr = 0.19;  end
if ~isfield(params,'Lt'), params.Lt = 0.514; end
if ~isfield(params,'Wt'), params.Wt = 0.639; end
if ~isfield(params,'m_r'), params.m_r = 50;  end
if ~isfield(params,'I_r'), params.I_r = 5;   end
if ~isfield(params,'m_t'), params.m_t = 20;  end
if ~isfield(params,'I_t'), params.I_t = 2;   end

% Initial state default 
if ~isfield(params,'x0')
    params.x0 = zeros(12,1);
end

% Project trailer initial state so the hitch points coincide exactly
xr0 = params.x0(1); yr0 = params.x0(2); thetar0 = params.x0(3);

% robot hitch in world
r_r_h = [xr0; yr0] - params.d*[cos(thetar0); sin(thetar0)];

% trailer COM must satisfy: r_th = r_r_h
% r_th = [xt; yt] + R(thetat) * [Lt/2 ; 0]
% choose initial thetat = thetar0
xt0 = r_r_h(1) - (params.Lt/2)*cos(thetar0);
yt0 = r_r_h(2) - (params.Lt/2)*sin(thetar0);
params.x0(7:9) = [xt0; yt0; thetar0];
params.x0(10:12) = [0;0;0];

% Integrate
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);
[t,X] = ode45(@(tt,xx) full_dynamics(tt,xx,controller_handle,params), ...
              tspan, params.x0, opts);

% Diagnostics
for k = 1:length(t)
    xr = X(k,1); yr = X(k,2); thetar = X(k,3);
    xt = X(k,7); yt = X(k,8); thetat = X(k,9);

    r_rh = [xr; yr] - params.d*[cos(thetar); sin(thetar)];
    r_th = [xt; yt] + (params.Lt/2)*[cos(thetat); sin(thetat)];

    hitch_len = norm(r_rh - r_th);
    fprintf('t=%.2f, hitch=%.4f\n', t(k), hitch_len);
end
% =========================================================
% 2) Build log: t, X, hitch output y, and control input u
% =========================================================
N = length(t);
log = struct();
log.t  = t;        % [N x 1] time vector
log.X  = X;        % [N x 12] full state trajectory
log.yd = zeros(N,2);
log.y  = zeros(N,2);   % actual hitch position (x, y)
log.u  = zeros(N,2);   % control inputs [F_drive, tau_r]
ref_state.mode = 'local_circle'; %% Notice: change this to draw different desire path
ref_state.initialized = false;
for k = 1:N
    % state and time at step k
    xk = X(k,:).';      % 12x1 state at time step k
    tk = t(k);

    % current hitch position (definition consistent with the controller)
    xr = xk(1); 
    yr = xk(2); 
    thetar = xk(3);
    d  = params.d;
    yk = [xr - d*cos(thetar);
          yr - d*sin(thetar)];
    log.y(k,:) = yk.';  % store actual hitch position

    % current control input (replay using the same controller_handle)
    [F_drive, tau_r] = controller_handle(xk, tk);
    log.u(k,:) = [F_drive, tau_r];

    % desired hitch trajectory evaluated at current time,
    % using actual hitch position as initial condition for hitch_ref
    [yd_k, ~, ~, ref_state] = hitch_ref(tk, yk, ref_state);
    log.yd(k,:) = yd_k.';    % store desired hitch position
end
save('log.mat','log')

% Visualization
figure('Color','w'); hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)');
title('Planar Robot–Trailer Dynamics (Rigid Hitch)');
xlim([-2 12]); ylim([-4 4]);

path_trace = plot(NaN,NaN,'k-','LineWidth',1.2,'HandleVisibility','off');
robot_patch   = patch(NaN,NaN,'r','FaceAlpha',0.4,'EdgeColor','none');
trailer_patch = patch(NaN,NaN,'b','FaceAlpha',0.4,'EdgeColor','none');

traj_x=[]; traj_y=[];
for i = 1:5:length(t)
    xr = X(i,1); yr = X(i,2); thetar = X(i,3);
    xt = X(i,7); yt = X(i,8); thetat = X(i,9);

    traj_x(end+1)=xr; traj_y(end+1)=yr;
    set(path_trace,'XData',traj_x,'YData',traj_y);

    set(robot_patch,'XData',rect_x(xr,params.Wr,params.Lr,thetar), ...
                    'YData',rect_y(yr,params.Wr,params.Lr,thetar));
    set(trailer_patch,'XData',rect_x(xt,params.Wt,params.Lt,thetat), ...
                      'YData',rect_y(yt,params.Wt,params.Lt,thetat));

    axis([xr-3 xr+5 yr-3 yr+3]);
    drawnow;
end
end

% Full rigid-body dynamics with exact constraints
function dx = full_dynamics(t, x, controller_handle, params)

xr = x(1); yr = x(2); thetar = x(3);
vxr = x(4); vyr = x(5); wr = x(6);
xt = x(7); yt = x(8); thetat = x(9);
vxt = x(10); vyt = x(11); wt = x(12);

m_r = params.m_r; I_r = params.I_r;
m_t = params.m_t; I_t = params.I_t;
d   = params.d;  Lt = params.Lt;

% Controller forces
[F_drive, tau_r] = controller_handle(x,t);
Fr_body = [cos(thetar); sin(thetar)] * F_drive;

% Hitch location offsets (body frame)
r_rh_body = [-d; 0];
r_th_body = [ Lt/2; 0];

% World hitch points
r_rh = [xr; yr] + R2(thetar) * r_rh_body;
r_th = [xt; yt] + R2(thetat) * r_th_body;

% Jacobians
p_r = R2(thetar) * r_rh_body;
p_t = R2(thetat) * r_th_body;
S_p_r = [-p_r(2); p_r(1)];
S_p_t = [-p_t(2); p_t(1)];

J_holo = [ eye(2),  S_p_r,  -eye(2), -S_p_t ];

l = [-sin(thetat); cos(thetat)];
J_nonholo = [0 0 0 l(1) l(2) 0];

J = [J_holo; J_nonholo];

M = diag([m_r m_r I_r m_t m_t I_t]);
v = [vxr; vyr; wr; vxt; vyt; wt];

Q = [Fr_body; tau_r; zeros(3,1)];

phi     = r_rh - r_th;
phi_dot = J_holo * v;

zeta = 0.9;
omega_b = 10.0;
b_holo = -(2*zeta*omega_b.*phi_dot + (omega_b^2).*phi);
b_nonholo = -J_nonholo * v;
b_total = [b_holo; b_nonholo];

A = [M, J'; J, zeros(size(J,1))];
rhs = [Q; b_total];
sol = A \ rhs;

a = sol(1:6);

fprintf('t=%.2f, hitch=%.4f\n', t, norm(phi));

dx = [vxr; vyr; wr; a(1); a(2); a(3); vxt; vyt; wt; a(4); a(5); a(6)];
end

% Helpers
function R = R2(th), R=[cos(th) -sin(th); sin(th) cos(th)]; end

function X = rect_x(xc,W,L,theta)
R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
pts=0.5*[-L -W; L -W; L W; -L W]';
X = R(1,:)*pts + xc;
end

function Y = rect_y(yc,W,L,theta)
R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
pts=0.5*[-L -W; L -W; L W; -L W]';
Y = R(2,:)*pts + yc;
end
