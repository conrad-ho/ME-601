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
[t,X] = ode45(@(tt,xx) ode_rhs_new_interface(tt, xx, controller_handle, params), ...
              tspan, params.x0, opts);

% diagnostics: hitch length check
for k = 1:length(t)
    xr = X(k,1); yr = X(k,2); thetar = X(k,3);
    xt = X(k,7); yt = X(k,8); thetat = X(k,9);

    r_rh = [xr; yr] - params.d*[cos(thetar); sin(thetar)];
    r_th = [xt; yt] + (params.Lt/2)*[cos(thetat); sin(thetat)];

    hitch_len = norm(r_rh - r_th);
    fprintf('t=%.2f, hitch=%.4f\n', t(k), hitch_len);
end

% =========================================================
% build log: t, X, trailer output y, and control input u
% note: here y follows the trailer com (x_t, y_t), not the hitch
% =========================================================
N = length(t);
log = struct();
log.t  = t;              % [N x 1]
log.X  = X;              % [N x 12]
log.yd = zeros(N,2);     % desired hitch position
log.y  = zeros(N,2);     % actual trailer com position
log.u  = zeros(N,2);     % control inputs [F_drive, tau_r]

ref_state.mode = 'local_circle';
ref_state.initialized = false;

for k = 1:N
    xk = X(k,:).';
    tk = t(k);

    xr     = xk(1);
    yr     = xk(2);
    thetar = xk(3);
    d      = params.d;

    y_hitch = [xr - d*cos(thetar);
               yr - d*sin(thetar)];

    xt = xk(7);
    yt = xk(8);
    y_trailer = [xt; yt];
    log.y(k,:) = y_trailer.';

    [F_drive, tau_r] = controller_handle(xk, tk);
    log.u(k,:) = [F_drive, tau_r];

    [yd_k, ~, ~, ref_state] = hitch_ref(tk, y_hitch, ref_state);
    log.yd(k,:) = yd_k.';
end

save('log.mat','log','-append')

% visualization
figure('Color','w'); hold on; grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)');
title('Planar Robot–Trailer Dynamics (Rigid Hitch)');
xlim([-2 12]); ylim([-4 4]);

path_robot   = plot(NaN,NaN,'k-','LineWidth',1.2,'DisplayName','robot');
path_trailer = plot(NaN,NaN,'k--','LineWidth',1.2,'DisplayName','trailer');

robot_patch   = patch(NaN,NaN,'r','FaceAlpha',0.4,'EdgeColor','none');
trailer_patch = patch(NaN,NaN,'b','FaceAlpha',0.4,'EdgeColor','none');

legend('Location','best');

idx_vec  = 1:5:length(t);
Nsteps   = numel(idx_vec);
traj_x_r = NaN(1, Nsteps); traj_y_r = NaN(1, Nsteps);
traj_x_t = NaN(1, Nsteps); traj_y_t = NaN(1, Nsteps);

kplot = 0;
for i = idx_vec
    kplot = kplot + 1;

    xr = X(i,1); yr = X(i,2); thetar = X(i,3);
    xt = X(i,7); yt = X(i,8); thetat = X(i,9);

    traj_x_r(kplot) = xr; traj_y_r(kplot) = yr;
    traj_x_t(kplot) = xt; traj_y_t(kplot) = yt;

    set(path_robot,  'XData', traj_x_r(1:kplot), 'YData', traj_y_r(1:kplot));
    set(path_trailer,'XData', traj_x_t(1:kplot), 'YData', traj_y_t(1:kplot));

    set(robot_patch,'XData',rect_x(xr,params.Wr,params.Lr,thetar), ...
                    'YData',rect_y(yr,params.Wr,params.Lr,thetar));
    set(trailer_patch,'XData',rect_x(xt,params.Wt,params.Lt,thetat), ...
                      'YData',rect_y(yt,params.Wt,params.Lt,thetat));

    axis([xr-3 xr+5 yr-3 yr+3]);
    drawnow;
end
end

% =========================================================
% ODE RHS using unified towing_dynamic interface
% =========================================================
function dx = ode_rhs_new_interface(t, x, controller_handle, params)
% x: 12x1, u: 2x1
% dx = towing_dynamic(x, params, 'full').xdot(u)

% controller input
[F_drive, tau_r] = controller_handle(x, t);
u = [F_drive; tau_r];

% unified dynamics
dyn = towing_dynamic(x, params, 'full');
dx  = dyn.xdot(u);
end

% helpers
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
