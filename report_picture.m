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
load('TO_Output.mat','to_log','to_dbg');
load('sim_log.mat','sim_log');
load('ref_to.mat','ref_to');
ref=ref_to;
% make sure ref hitch is Nx2
    p_h_ref = ref.p_h;
    if size(p_h_ref,1) == 2
        p_h_ref = p_h_ref.';
    end


% ---------------- figure 1: hitch xy (ref vs TO vs QP execution) ----------------
figure('Color','w'); hold on; grid on; axis equal;
title('hitch trajectory: reference vs TO plan vs QP execution');
xlabel('x (m)'); ylabel('y (m)');

% 1) reference hitch
plot(p_h_ref(:,1), p_h_ref(:,2), 'k:', 'LineWidth', 1.8);

% 2) TO plan hitch (from to_log)
if isfield(to_log,'p_h_r') && ~isempty(to_log.p_h_r)
    plot(to_log.p_h_r(:,1), to_log.p_h_r(:,2), 'b--', 'LineWidth', 1.6);
end

% 3) QP execution hitch (from sim_log)
if isfield(sim_log,'p_h_r') && ~isempty(sim_log.p_h_r)
    plot(sim_log.p_h_r(:,1), sim_log.p_h_r(:,2), 'r-', 'LineWidth', 1.8);
end

legend({'reference (hitch)', 'TO plan (hitch)', 'QP execution (hitch)'}, ...
       'Location', 'best');
% --- make figure size for ieee single-column ---
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 3.5 2.6]);  % 宽 3.5in 单栏常用；高按内容调

% 字体别太大也别太小
set(gca, 'FontSize', 9);
exportgraphics(gcf, fullfile('figures','fig_hitch_xy.pdf'), ...
    'ContentType','vector', 'Resolution', 300);
%%
% ---------------- figure 2: hitch tracking error vs time ----------------
figure('Color','w'); hold on; grid on;
title('hitch tracking error vs time');
xlabel('time (s)'); ylabel('||e|| (m)');

% --- time vectors ---
t_sim = sim_log.t(:);          % 3001x1
t_ref = ref_to.t(:);           % 151x1
t_to  = to_log.t(:);           % 151x1

% --- hitch positions (make them Nx2) ---
% ref_to.p_h is 2x151 (from your screenshot)
p_ref = ref_to.p_h;
if size(p_ref,1) == 2
    p_ref = p_ref.';           % -> 151x2
end

% to_log planned hitch (151x2)
p_to = to_log.p_h_r;

% sim executed hitch (3001x2)
p_sim = sim_log.p_h_r;

% --- interpolate ref and TO to sim time grid (assume same start/end time) ---
p_ref_i = interp1(t_ref, p_ref, t_sim, 'linear', 'extrap');
p_to_i  = interp1(t_to,  p_to,  t_sim, 'linear', 'extrap');

% --- errors ---
e_qp = p_sim   - p_ref_i;      % QP execution vs reference
e_to = p_to_i  - p_ref_i;      % TO plan vs reference

err_qp = vecnorm(e_qp, 2, 2);  % 3001x1
err_to = vecnorm(e_to, 2, 2);  % 3001x1

% --- plot (keep it clean) ---
plot(t_sim, err_to, 'b--', 'LineWidth', 1.4);
plot(t_sim, err_qp, 'r-',  'LineWidth', 1.6);

legend({'TO plan error ||y_{to}-y_{ref}||', ...
        'QP execution error ||y_{sim}-y_{ref}||'}, 'Location', 'best');

% --- optional: annotate RMS / max (uses to_dbg if you want) ---
if exist('to_dbg','var') && ~isempty(to_dbg)
    if isfield(to_dbg,'hitch_err_rms') && isfield(to_dbg,'hitch_err_max')
        txt = sprintf('TO rms=%.3f m, max=%.3f m', to_dbg.hitch_err_rms, to_dbg.hitch_err_max);
        x0 = t_sim( round(0.60*numel(t_sim)) );
        y0 = max([err_qp; err_to]) * 0.90;
        text(x0, y0, txt, 'FontSize', 9);
    end
end

% --- make figure size for ieee single-column + export ---
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 3.5 2.6]);  % single-column friendly
set(gca, 'FontSize', 9);

% export (vector pdf)
if ~exist('figures','dir'), mkdir('figures'); end
exportgraphics(gcf, fullfile('figures','fig2_hitch_error.pdf'), ...
    'ContentType','vector', 'Resolution', 300);
%%
% ---------------- figure 3: Lyapunov function V(t) ----------------
figure('Color','w'); hold on; grid on;
title('Lyapunov function for hitch tracking');
xlabel('time (s)'); ylabel('V(t)');

t_sim = sim_log.t(:);
p_sim = sim_log.p_h_r;   % 3001x2

% ref hitch position (2x151 -> 151x2)
p_ref = ref_to.p_h;
if size(p_ref,1) == 2
    p_ref = p_ref.';
end

% ref hitch velocity (2x151 -> 151x2)
v_ref = ref_to.pdot_h;
if size(v_ref,1) == 2
    v_ref = v_ref.';
end

t_ref = ref_to.t(:);

% interpolate reference to sim time grid
p_ref_i = interp1(t_ref, p_ref, t_sim, 'linear', 'extrap');
v_ref_i = interp1(t_ref, v_ref, t_sim, 'linear', 'extrap');

% compute sim hitch velocity (numerical derivative)
dt = mean(diff(t_sim));
v_sim = gradient(p_sim, dt);     % 3001x2 (more stable than diff)

% errors
e  = p_sim - p_ref_i;            % 3001x2
ed = v_sim - v_ref_i;            % 3001x2

% pick Kp (use your actual Kp if available)
if isfield(params,'Kp_tr')
    Kp = params.Kp_tr;
else
    Kp = diag([10, 10]);         % fallback (must be PD)
end

% Lyapunov V(t) = 1/2||ed||^2 + 1/2 e'Kp e
V = 0.5*sum(ed.^2,2) + 0.5*sum((e*Kp).*e, 2);

plot(t_sim, V, 'k-', 'LineWidth', 1.6);

% optional: show Vdot (numerical) as dashed line if you want
% Vdot = gradient(V, dt);
% plot(t_sim, Vdot, 'k--', 'LineWidth', 1.0);

% ieee single-column sizing + export
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [1 1 3.5 2.6]);
set(gca, 'FontSize', 9);

if ~exist('figures','dir'), mkdir('figures'); end
exportgraphics(gcf, fullfile('figures','fig3_lyapunov.pdf'), ...
    'ContentType','vector', 'Resolution', 300);
