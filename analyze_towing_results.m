function analyze_towing_results(sim, to, ref, params)
% analyze_towing_results (ref is hitch-only)
% plots:
%   1) xy path compare (robot/trailer + hitch, ref hitch)
%   2) yaw compare (sim vs to; optional ref theta_h)
%   3) hitch overlap (constraint) + hitch tracking (vs ref)

    % make sure ref hitch is Nx2
    p_h_ref = ref.p_h;
    if size(p_h_ref,1) == 2
        p_h_ref = p_h_ref.';
    end

    % ---------------- figure 1: xy path compare ----------------
    figure('Color','w'); hold on; grid on; axis equal;
    title('path compare: sim vs to vs ref (ref = hitch only)');
    xlabel('x (m)'); ylabel('y (m)');

    % ref hitch
    plot(p_h_ref(:,1), p_h_ref(:,2), 'k:', 'LineWidth',1.6);

    % to
    plot(to.p_r(:,1),  to.p_r(:,2),  'b-',  'LineWidth',1.2);
    plot(to.p_tr(:,1), to.p_tr(:,2), 'b--', 'LineWidth',1.2);
    if isfield(to,'p_h_r') && ~isempty(to.p_h_r)
        plot(to.p_h_r(:,1), to.p_h_r(:,2), 'b:', 'LineWidth',1.1);
    end

    % sim
    plot(sim.p_r(:,1),  sim.p_r(:,2),  'r-',  'LineWidth',1.4);
    plot(sim.p_tr(:,1), sim.p_tr(:,2), 'r--', 'LineWidth',1.4);
    if isfield(sim,'p_h_r') && ~isempty(sim.p_h_r)
        plot(sim.p_h_r(:,1), sim.p_h_r(:,2), 'r:', 'LineWidth',1.2);
    end

    legend({'ref hitch','to robot','to trailer','to hitch', ...
            'sim robot','sim trailer','sim hitch'}, 'Location','best');

    % ---------------- figure 2: yaw compare ----------------
    figure('Color','w'); grid on; hold on;
    title('yaw compare (robot yaw)');
    xlabel('t (s)'); ylabel('theta (rad)');

    plot(to.t,  unwrap(to.theta),  'b-', 'LineWidth',1.2);
    plot(sim.t, unwrap(sim.theta), 'r-', 'LineWidth',1.4);

    % optional: show hitch tangent heading if provided
    if isfield(ref,'theta_h') && ~isempty(ref.theta_h)
        plot(ref.t(:), unwrap(ref.theta_h(:)), 'k:', 'LineWidth',1.2);
        legend({'to theta_r','sim theta_r','ref theta_h'}, 'Location','best');
    else
        legend({'to theta_r','sim theta_r'}, 'Location','best');
    end

    % integrated yaw from wr (helps diagnose drift)
    if isfield(sim,'wr') && ~isempty(sim.wr)
        theta_int_sim = sim.theta(1) + cumtrapz(sim.t, sim.wr);
        plot(sim.t, unwrap(theta_int_sim), 'r--', 'LineWidth',1.0);
    end
    if isfield(to,'wr') && ~isempty(to.wr)
        theta_int_to = to.theta(1) + cumtrapz(to.t, to.wr);
        plot(to.t, unwrap(theta_int_to), 'b--', 'LineWidth',1.0);
    end

    % ---------------- figure 3: hitch metrics ----------------
    figure('Color','w'); grid on; hold on;
    title('hitch metrics');
    xlabel('t (s)'); ylabel('meters');

    % (a) constraint overlap error: ||p_h_r - p_h_t||
    plot(to.t,  to.e_h,  'b-', 'LineWidth',1.2);
    plot(sim.t, sim.e_h, 'r-', 'LineWidth',1.4);

    % (b) hitch tracking error vs ref: ||p_h_r - p_h_ref||
    % interpolate ref hitch onto sim/to time grids
    p_h_ref_to  = interp1(ref.t(:), p_h_ref, to.t(:),  'linear', 'extrap');
    p_h_ref_sim = interp1(ref.t(:), p_h_ref, sim.t(:), 'linear', 'extrap');

    e_track_to  = sqrt(sum((to.p_h_r  - p_h_ref_to ).^2, 2));
    e_track_sim = sqrt(sum((sim.p_h_r - p_h_ref_sim).^2, 2));

    plot(to.t,  e_track_to,  'b--', 'LineWidth',1.0);
    plot(sim.t, e_track_sim, 'r--', 'LineWidth',1.0);

    legend({'to hitch overlap (constraint)', 'sim hitch overlap (constraint)', ...
            'to hitch tracking (vs ref)',     'sim hitch tracking (vs ref)'}, ...
            'Location','best');

    % ---------------- optional: rv if available ----------------
    if isfield(sim,'rv') && ~isempty(sim.rv)
        figure('Color','w'); grid on; hold on;
        title('constraint velocity residual  ||Jv|| (sim)');
        xlabel('t (s)'); ylabel('rv');
        plot(sim.t, sim.rv, 'r-', 'LineWidth',1.3);
    end
end
