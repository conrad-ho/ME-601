function check_to_constraints(x_to, params, dt_to)
%CHECK_TO_CONSTRAINTS
%   Check constraint quality along a TO trajectory:
%   - hitch position error |Δp_h|
%   - constraint velocity residual ||J v||
%   - Baumgarte residual ||J a + dotJ*v + 2α J v + β^2 h||
%
% x_to   : [12 x (N+1)] TO state trajectory
% params : struct (must contain d, Lt, and optionally alpha_baum, beta_baum)
% dt_to  : nominal time step between nodes

    if ~isfield(params, 'alpha_baum'), params.alpha_baum = 5.0;  end
    if ~isfield(params, 'beta_baum'),  params.beta_baum  = 10.0; end
    alpha = params.alpha_baum;
    beta  = params.beta_baum;

    d  = params.d;
    Lt = params.Lt;
    x_to=x_to.';
    [nx, Np1] = size(x_to);
    if nx ~= 12
        error('x_to must be 12 x (N+1).');
    end
    N = Np1 - 1;

    phi_norm   = zeros(Np1,1);  % hitch position error
    Jv_norm    = zeros(Np1,1);  % ||J v||
    r_baum_norm = zeros(N,1);   % Baumgarte residual

    % === loop over all nodes to compute |Δp_h| and ||J v|| ===
    for k = 1:Np1
        xk = x_to(:,k);

        xr    = xk(1);  yr    = xk(2);  thetar  = xk(3);
        vxr   = xk(4);  vyr   = xk(5);  wr      = xk(6);
        xt    = xk(7);  yt    = xk(8);  thetat  = xk(9);
        vxt   = xk(10); vyt   = xk(11); wt      = xk(12);

        % hitch offsets in body frame
        r_rh_b = [-d;   0];
        r_th_b = [ Lt/2; 0];

        Rr = R2(thetar);
        Rt = R2(thetat);

        p_r = Rr * r_rh_b;
        p_t = Rt * r_th_b;

        p_h_r  = [xr; yr] + p_r;
        p_h_tr = [xt; yt] + p_t;

        phi = p_h_r - p_h_tr;          % 2x1
        phi_norm(k) = norm(phi);

        % J*v
        dyn = towing_dynamics_mats(xk, params);
        J   = dyn.J;                   % 3x6
        v   = dyn.v;                   % 6x1
        Jv  = J * v;
        Jv_norm(k) = norm(Jv);
    end

    % === approximate a by finite differences and check Baumgarte residual ===
    for k = 1:N
        xk   = x_to(:,k);
        xkp1 = x_to(:,k+1);

        % velocities at k and k+1
        vk   = [xk(4:6);   xk(10:12)];    % 6x1
        vk1  = [xkp1(4:6); xkp1(10:12)];  % 6x1

        a_fd = (vk1 - vk) / dt_to;        % finite-diff approx of a

        xr    = xk(1);  yr    = xk(2);  thetar  = xk(3);
        xt    = xk(7);  yt    = xk(8);  thetat  = xk(9);

        % constraint again to build h
        r_rh_b = [-d;   0];
        r_th_b = [ Lt/2; 0];

        Rr = R2(thetar);
        Rt = R2(thetat);

        p_r = Rr * r_rh_b;
        p_t = Rt * r_th_b;

        p_h_r  = [xr; yr] + p_r;
        p_h_tr = [xt; yt] + p_t;

        phi = p_h_r - p_h_tr;          % 2x1
        h   = [phi;
               0];                     % third constraint has no position level

        % dynamics matrices
        dyn = towing_dynamics_mats(xk, params);
        J   = dyn.J;
        dotJ= dyn.dotJ;
        v   = dyn.v;

        % Baumgarte equation residual:
        % r = J*a + dotJ*v + 2 alpha J v + beta^2 h
        Jv  = J * v;
        r   = J * a_fd + dotJ * v + 2*alpha*Jv + (beta^2)*h;

        r_baum_norm(k) = norm(r);
    end

    % === time vector ===
    t = (0:Np1-1).' * dt_to;
    t_mid = (0:N-1).' * dt_to + dt_to/2;  % for residuals between nodes

    % === plots ===
    figure;
    subplot(3,1,1);
    plot(t, phi_norm, 'LineWidth', 1.5);
    ylabel('|Δp_h| (m)');
    grid on;
    title('TO trajectory constraint diagnostics');

    subplot(3,1,2);
    plot(t, Jv_norm, 'LineWidth', 1.5);
    ylabel('||J v||');
    grid on;

    subplot(3,1,3);
    plot(t_mid, r_baum_norm, 'LineWidth', 1.5);
    ylabel('||Baumgarte residual||');
    xlabel('time (s)');
    grid on;
end
function R = R2(theta)
%R2   2×2 rotation matrix for angle theta (radians)
%
%   R = [ cos(theta), -sin(theta);
%         sin(theta),  cos(theta) ];

    R = [cos(theta), -sin(theta);
         sin(theta),  cos(theta)];
end