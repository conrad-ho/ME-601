function test_normalized_dynamics_constraints()
%TEST_NORMALIZED_DYNAMICS_CONSTRAINTS
%   Sanity check for the normalized towing_dynamics_full and
%   towing_dynamics_mats.
%
%   - Builds a constraint-consistent initial state x0
%   - Simulates with a simple controller
%   - Plots |Δp_h| and ||J v|| over time

    % === parameters ===
    params.m_r = 50;
    params.I_r = 2;
    params.m_t = 20;
    params.I_t = 1;
    params.d   = 0.134;
    params.Lt  = 0.514;

    % Baumgarte gains
    params.alpha_baum = 5.0;
    params.beta_baum  = 10.0;

    % === build constraint-consistent x0 ===
    xr0     = 0.0;
    yr0     = 0.0;
    thetar0 = 0.0;

    thetat0 = 0.0;

    r_rh_b = [-params.d;   0];
    r_th_b = [ params.Lt/2; 0];

    Rr0 = [cos(thetar0), -sin(thetar0);
           sin(thetar0),  cos(thetar0)];
    Rt0 = [cos(thetat0), -sin(thetat0);
           sin(thetat0),  cos(thetat0)];

    p_r0 = Rr0 * r_rh_b;
    p_t0 = Rt0 * r_th_b;

    xt0 = xr0 + p_r0(1) - p_t0(1);
    yt0 = yr0 + p_r0(2) - p_t0(2);

    vxr0 = 0;  vyr0 = 0;  wr0 = 0;
    vxt0 = 0;  vyt0 = 0;  wt0 = 0;

    x0 = [xr0; yr0; thetar0; vxr0; vyr0; wr0;
          xt0; yt0; thetat0; vxt0; vyt0; wt0];

    % === simple controller (here zero input) ===
    ctrl = @(x,t) [0; 0];   % you can change to any (F_drive,tau_r)

    % === simulate ===
    tspan = [0 10];
    odefun = @(t,x) towing_dynamics_full(x, ctrl(x,t), params);
    opts   = odeset('RelTol',1e-8, 'AbsTol',1e-10);

    [t_sol, X_sol] = ode45(odefun, tspan, x0, opts);

    % === post-process: |Δp_h| and ||J v|| ===
    n = numel(t_sol);
    phi_norm = zeros(n,1);
    Jv_norm  = zeros(n,1);

    for k = 1:n
        xk = X_sol(k,:).';

        xr = xk(1);  yr = xk(2);  thetar = xk(3);
        xt = xk(7);  yt = xk(8);  thetat = xk(9);

        % hitch positions (same as in full dynamics)
        r_rh_b = [-params.d;   0];
        r_th_b = [ params.Lt/2; 0];

        Rr = [cos(thetar), -sin(thetar);
              sin(thetar),  cos(thetar)];
        Rt = [cos(thetat), -sin(thetat);
              sin(thetat),  cos(thetat)];

        p_r = Rr * r_rh_b;
        p_t = Rt * r_th_b;

        p_h_r  = [xr; yr] + p_r;
        p_h_tr = [xt; yt] + p_t;

        phi = p_h_r - p_h_tr;
        phi_norm(k) = norm(phi);

        % J v
        dyn = towing_dynamics_mats(xk, params);
        J   = dyn.J;
        v   = dyn.v;
        Jv  = J * v;
        Jv_norm(k) = norm(Jv);
    end

    % === plots ===
    figure;
    subplot(2,1,1);
    plot(t_sol, phi_norm, 'LineWidth', 1.5);
    ylabel('|Δp_h| (m)');
    grid on;
    title('Constraint diagnostics for normalized dynamics');

    subplot(2,1,2);
    plot(t_sol, Jv_norm, 'LineWidth', 1.5);
    ylabel('||J v||');
    xlabel('time (s)');
    grid on;
end
