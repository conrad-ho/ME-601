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
        % ---------------- figure 4: qp hitch monitors (h/hdot/hddot) ----------------
    % expects (optional) fields:
    %   sim.qp_h [N x 2], sim.qp_hdot [N x 2], sim.qp_e_yddot [N x 2] (or sim.qp_yddot)
    %   to.qp_h  [N x 2], to.qp_hdot  [N x 2], to.qp_e_yddot  [N x 2] (or to.qp_yddot)

    has_sim_h    = isfield(sim,'qp_h')       && ~isempty(sim.qp_h);
    has_sim_hdot = isfield(sim,'qp_hdot')    && ~isempty(sim.qp_hdot);
    has_sim_ea   = isfield(sim,'qp_e_yddot') && ~isempty(sim.qp_e_yddot);
    has_sim_ya   = isfield(sim,'qp_yddot')   && ~isempty(sim.qp_yddot);

    has_to_h     = isfield(to,'qp_h')        && ~isempty(to.qp_h);
    has_to_hdot  = isfield(to,'qp_hdot')     && ~isempty(to.qp_hdot);
    has_to_ea    = isfield(to,'qp_e_yddot')  && ~isempty(to.qp_e_yddot);
    has_to_ya    = isfield(to,'qp_yddot')    && ~isempty(to.qp_yddot);

    if has_sim_h || has_sim_hdot || has_sim_ea || has_sim_ya || has_to_h || has_to_hdot || has_to_ea || has_to_ya
        figure('Color','w');
        tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

        % --- (1) ||h|| ---
        nexttile; grid on; hold on;
        title('qp monitor: ||h|| (h = y - y_d)');
        xlabel('t (s)'); ylabel('||h||');

        if has_to_h
            plot(to.t, sqrt(sum(to.qp_h.^2,2)), 'b-', 'LineWidth',1.1);
        end
        if has_sim_h
            plot(sim.t, sqrt(sum(sim.qp_h.^2,2)), 'r-', 'LineWidth',1.3);
        end
        legend_entries = {};
        if has_to_h,  legend_entries{end+1}  = 'to ||h||';  end
        if has_sim_h, legend_entries{end+1}  = 'sim ||h||'; end
        if ~isempty(legend_entries), legend(legend_entries,'Location','best'); end

        % --- (2) ||hdot|| ---
        nexttile; grid on; hold on;
        title('qp monitor: ||hdot|| (hdot = ydot - ydot_d)');
        xlabel('t (s)'); ylabel('||hdot||');

        if has_to_hdot
            plot(to.t, sqrt(sum(to.qp_hdot.^2,2)), 'b-', 'LineWidth',1.1);
        end
        if has_sim_hdot
            plot(sim.t, sqrt(sum(sim.qp_hdot.^2,2)), 'r-', 'LineWidth',1.3);
        end
        legend_entries = {};
        if has_to_hdot,  legend_entries{end+1}  = 'to ||hdot||';  end
        if has_sim_hdot, legend_entries{end+1}  = 'sim ||hdot||'; end
        if ~isempty(legend_entries), legend(legend_entries,'Location','best'); end

        % --- (3) ||e_yddot|| preferred, else ||yddot|| ---
        nexttile; grid on; hold on;
        if (has_to_ea || has_sim_ea)
            title('qp monitor: ||e_{yddot}|| (e = yddot - yddot_{des})');
            ylabel('||e_{yddot}||');
            if has_to_ea
                plot(to.t, sqrt(sum(to.qp_e_yddot.^2,2)), 'b-', 'LineWidth',1.1);
            end
            if has_sim_ea
                plot(sim.t, sqrt(sum(sim.qp_e_yddot.^2,2)), 'r-', 'LineWidth',1.3);
            end
            legend_entries = {};
            if has_to_ea,  legend_entries{end+1}  = 'to ||e_{yddot}||';  end
            if has_sim_ea, legend_entries{end+1}  = 'sim ||e_{yddot}||'; end
        else
            title('qp monitor: ||yddot|| (actual)');
            ylabel('||yddot||');
            if has_to_ya
                plot(to.t, sqrt(sum(to.qp_yddot.^2,2)), 'b-', 'LineWidth',1.1);
            end
            if has_sim_ya
                plot(sim.t, sqrt(sum(sim.qp_yddot.^2,2)), 'r-', 'LineWidth',1.3);
            end
            legend_entries = {};
            if has_to_ya,  legend_entries{end+1}  = 'to ||yddot||';  end
            if has_sim_ya, legend_entries{end+1}  = 'sim ||yddot||'; end
        end
        xlabel('t (s)');
        if ~isempty(legend_entries), legend(legend_entries,'Location','best'); end
    end
    
    % ---------------- optional: rv if available ----------------
    if isfield(sim,'rv') && ~isempty(sim.rv)
        figure('Color','w'); grid on; hold on;
        title('constraint velocity residual  ||Jv|| (sim)');
        xlabel('t (s)'); ylabel('rv');
        plot(sim.t, sim.rv, 'r-', 'LineWidth',1.3);
    end
end
