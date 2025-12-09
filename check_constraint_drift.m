function check_constraint_drift()
%CHECK_CONSTRAINT_DRIFT
%   Sanity test for Baumgarte-stabilized rigid hitch dynamics.
%   - builds a constraint-consistent initial state x0
%   - simulates full_dynamics for some time
%   - plots |delta p_h| and |J v| over time

    % ==== 1) set parameters ====
    params.m_r = 50;    % robot mass
    params.I_r = 2;     % robot inertia
    params.m_t = 20;    % trailer mass
    params.I_t = 1;     % trailer inertia
    params.d   = 0.134; % hitch offset on robot
    params.Lt  = 0.514; % trailer length (or COM offset definition)

    % Baumgarte gains
    params.alpha_baum = 5.0;
    params.beta_baum  = 10.0;

    % ==== 2) build a constraint-consistent initial state x0 ====
    % We choose a simple configuration where robot and trailer are aligned.

    xr0     = 0.0;
    yr0     = 0.0;
    thetar0 = 0.0;   % robot facing +x

    thetat0 = 0.0;   % trailer aligned with robot

    % hitch offsets in body frames
    r_rh_b = [-params.d;   0];
    r_th_b = [ params.Lt/2; 0];

    Rr0    = R2(thetar0);
    Rt0    = R2(thetat0);

    p_r0   = Rr0 * r_rh_b;
    p_t0   = Rt0 * r_th_b;

    % enforce hitch coincidence: xr + p_r = xt + p_t
    xt0 = xr0 + p_r0(1) - p_t0(1);
    yt0 = yr0 + p_r0(2) - p_t0(2);

    % zero velocities (this automatically satisfies J*v = 0)
    vxr0 = 0;  vyr0 = 0;  wr0 = 0;
    vxt0 = 0;  vyt0 = 0;  wt0 = 0;

    x0 = [xr0; yr0; thetar0; vxr0; vyr0; wr0;
          xt0; yt0; thetat0; vxt0; vyt0; wt0];

    % ==== 3) choose a simple controller for the test ====
    % You can switch between:
    %   - zero input:         no_drive_controller
    %   - simple constant u:  simple_drive_controller

    controller_handle = @(x,t) no_drive_controller(x,t);
    % controller_handle = @(x,t) simple_drive_controller(x,t);

    % ==== 4) simulate with full_dynamics ====
    tspan = [0 10];  % simulate 10 seconds

    odefun = @(t,x) full_dynamics(t, x, controller_handle, params);
    opts   = odeset('RelTol',1e-8, 'AbsTol',1e-10);

    [t_sol, X_sol] = ode45(odefun, tspan, x0, opts);

    % ==== 5) post-process: compute |delta p_h| and |J v| ====
    n = numel(t_sol);
    phi_norm = zeros(n,1);
    Jv_norm  = zeros(n,1);

    for k = 1:n
        xk = X_sol(k,:).';

        xr    = xk(1);  yr    = xk(2);  thetar  = xk(3);
        vxr   = xk(4);  vyr   = xk(5);  wr      = xk(6);
        xt    = xk(7);  yt    = xk(8);  thetat  = xk(9);
        vxt   = xk(10); vyt   = xk(11); wt      = xk(12);

        % hitch offsets
        r_rh_b = [-params.d;   0];
        r_th_b = [ params.Lt/2; 0];

        Rr = R2(thetar);
        Rt = R2(thetat);

        p_r = Rr * r_rh_b;
        p_t = Rt * r_th_b;

        p_h_r  = [xr; yr] + p_r;
        p_h_tr = [xt; yt] + p_t;

        phi = p_h_r - p_h_tr;   % holonomic position error
        phi_norm(k) = norm(phi);

        % J*v
        dyn = towing_dynamics_mats(xk, params);
        J   = dyn.J;
        v   = dyn.v;
        Jv  = J * v;
        Jv_norm(k) = norm(Jv);
    end

    % ==== 6) plot results ====
    figure;
    subplot(2,1,1);
    plot(t_sol, phi_norm, 'LineWidth', 1.5);
    ylabel('|Δp_h| (m)');
    grid on;
    title('Hitch position error and constraint velocity residual');

    subplot(2,1,2);
    plot(t_sol, Jv_norm, 'LineWidth', 1.5);
    ylabel('||J v||');
    xlabel('time (s)');
    grid on;

end

% === helper controllers for testing ===

function [F_drive, tau_r] = no_drive_controller(~, ~)
%NO_DRIVE_CONTROLLER  zero input (pure drift test)
    F_drive = 0.0;
    tau_r   = 0.0;
end

function [F_drive, tau_r] = simple_drive_controller(~, t)
%SIMPLE_DRIVE_CONTROLLER  small drive + yaw pulse
    F_drive = 20.0;
    if t < 1.0
        tau_r = 0.0;
    else
        tau_r = 0.5;
    end
end
function R = R2(theta)
%R2   2×2 rotation matrix for angle theta (radians)
%
%   R = [ cos(theta), -sin(theta);
%         sin(theta),  cos(theta) ];

    R = [cos(theta), -sin(theta);
         sin(theta),  cos(theta)];
end
function dx = full_dynamics(t, x, controller_handle, params)
% full rigid-body dynamics with baumgarte-stabilized constraints
%
% x = [xr; yr; thetar; vxr; vyr; wr; xt; yt; thetat; vxt; vyt; wt]

    % unpack state (for clarity)
    xr     = x(1);  yr     = x(2);  thetar  = x(3);
    vxr    = x(4);  vyr    = x(5);  wr      = x(6);
    xt     = x(7);  yt     = x(8);  thetat  = x(9);
    vxt    = x(10); vyt    = x(11); wt      = x(12);

    % baumgarte gains
    if ~isfield(params, 'alpha_baum'), params.alpha_baum = 5.0;  end
    if ~isfield(params, 'beta_baum'),  params.beta_baum  = 10.0; end
    alpha = params.alpha_baum;
    beta  = params.beta_baum;

    d  = params.d;
    Lt = params.Lt;

    % controller inputs
    [F_drive, tau_r] = controller_handle(x, t);

    % generalized dynamics matrices (reuse same helper as TO)
    dyn  = towing_dynamics_mats(x, params);
    M    = dyn.M;       % 6x6
    B    = dyn.B;       % 6x2
    J    = dyn.J;       % 3x6
    dotJ = dyn.dotJ;    % 3x6
    v    = dyn.v;       % 6x1

    % generalized input forces
    u   = [F_drive; tau_r];
    Q   = B * u;        % 6x1

    % hitch position error h(q) = [phi; 0]
    % world hitch offsets (same geometry as in towing_dynamics_mats)
    r_rh_body = [-d; 0];
    r_th_body = [ Lt/2; 0];

    Rr = R2(thetar);
    Rt = R2(thetat);

    p_r = Rr * r_rh_body;           % hitch vector from robot com
    p_t = Rt * r_th_body;           % hitch vector from trailer com

    p_h_r  = [xr; yr] + p_r;        % hitch on robot
    p_h_tr = [xt; yt] + p_t;        % hitch on trailer

    phi = p_h_r - p_h_tr;           % 2x1 holonomic position error

    % third constraint (nonholonomic) has no position level -> set to 0
    h = [phi;
         0];                        % 3x1

    % baumgarte-stabilized constraint rhs:
    % J a + dotJ v + 2 alpha J v + beta^2 h = 0
    Jv      = J * v;
    rhs_con = -dotJ * v - 2*alpha*Jv - (beta^2)*h;   % 3x1

    % kkt system:
    % [M  J';  J  0] [a; lambda] = [Q; rhs_con]
    A   = [M,  J.';
           J,  zeros(size(J,1))];
    rhs = [Q;
           rhs_con];

    sol = A \ rhs;
    a   = sol(1:6);

    % assemble state derivative
    dx = zeros(12,1);
    dx(1)  = vxr;
    dx(2)  = vyr;
    dx(3)  = wr;
    dx(4)  = a(1);
    dx(5)  = a(2);
    dx(6)  = a(3);
    dx(7)  = vxt;
    dx(8)  = vyt;
    dx(9)  = wt;
    dx(10) = a(4);
    dx(11) = a(5);
    dx(12) = a(6);
end