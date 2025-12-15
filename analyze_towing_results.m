function analyze_towing_results(sim_log, to_log, ref, params, to_dbg)
% analyze_towing_results (ref is hitch-only)
% plots:
%   1) xy path compare (robot/trailer + hitch, ref hitch)
%   2) yaw compare (sim vs to; optional ref theta_h)
%   3) hitch overlap (constraint) + hitch tracking (vs ref)
%   4) qp hitch monitors (optional)
%   5) input margin (optional)
%   6) to cost breakdown (optional, needs to_dbg)

    if nargin < 5
        to_dbg = [];
    end

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

    % to (to_log)
    if isfield(to_log,'p_r') && ~isempty(to_log.p_r)
        plot(to_log.p_r(:,1),  to_log.p_r(:,2),  'b-',  'LineWidth',1.2);
    end
    if isfield(to_log,'p_tr') && ~isempty(to_log.p_tr)
        plot(to_log.p_tr(:,1), to_log.p_tr(:,2), 'b--', 'LineWidth',1.2);
    end
    if isfield(to_log,'p_h_r') && ~isempty(to_log.p_h_r)
        plot(to_log.p_h_r(:,1), to_log.p_h_r(:,2), 'b:', 'LineWidth',1.1);
    end

    % sim (sim_log)
    if isfield(sim_log,'p_r') && ~isempty(sim_log.p_r)
        plot(sim_log.p_r(:,1),  sim_log.p_r(:,2),  'r-',  'LineWidth',1.4);
    end
    if isfield(sim_log,'p_tr') && ~isempty(sim_log.p_tr)
        plot(sim_log.p_tr(:,1), sim_log.p_tr(:,2), 'r--', 'LineWidth',1.4);
    end
    if isfield(sim_log,'p_h_r') && ~isempty(sim_log.p_h_r)
        plot(sim_log.p_h_r(:,1), sim_log.p_h_r(:,2), 'r:', 'LineWidth',1.2);
    end

    legend({'ref hitch','to robot','to trailer','to hitch', ...
            'sim robot','sim trailer','sim hitch'}, 'Location','best');

    % ---------------- figure 2: yaw compare ----------------
    figure('Color','w'); grid on; hold on;
    title('yaw compare (robot yaw)');
    xlabel('t (s)'); ylabel('theta (rad)');

    if isfield(to_log,'t') && isfield(to_log,'theta') && ~isempty(to_log.t) && ~isempty(to_log.theta)
        plot(to_log.t,  unwrap(to_log.theta),  'b-', 'LineWidth',1.2);
    end
    if isfield(sim_log,'t') && isfield(sim_log,'theta') && ~isempty(sim_log.t) && ~isempty(sim_log.theta)
        plot(sim_log.t, unwrap(sim_log.theta), 'r-', 'LineWidth',1.4);
    end

    % optional: show hitch tangent heading if provided
    if isfield(ref,'theta_h') && ~isempty(ref.theta_h) && isfield(ref,'t') && ~isempty(ref.t)
        plot(ref.t(:), unwrap(ref.theta_h(:)), 'k:', 'LineWidth',1.2);
        legend({'to theta_r','sim theta_r','ref theta_h'}, 'Location','best');
    else
        legend({'to theta_r','sim theta_r'}, 'Location','best');
    end

    % integrated yaw from wr (helps diagnose drift)
    if isfield(sim_log,'wr') && ~isempty(sim_log.wr) && isfield(sim_log,'t') && ~isempty(sim_log.t) ...
            && isfield(sim_log,'theta') && ~isempty(sim_log.theta)
        theta_int_sim = sim_log.theta(1) + cumtrapz(sim_log.t, sim_log.wr);
        plot(sim_log.t, unwrap(theta_int_sim), 'r--', 'LineWidth',1.0);
    end
    if isfield(to_log,'wr') && ~isempty(to_log.wr) && isfield(to_log,'t') && ~isempty(to_log.t) ...
            && isfield(to_log,'theta') && ~isempty(to_log.theta)
        theta_int_to = to_log.theta(1) + cumtrapz(to_log.t, to_log.wr);
        plot(to_log.t, unwrap(theta_int_to), 'b--', 'LineWidth',1.0);
    end

    % ---------------- figure 3: hitch metrics ----------------
    figure('Color','w'); grid on; hold on;
    title('hitch metrics');
    xlabel('t (s)'); ylabel('meters');

    % (a) constraint overlap error: ||p_h_r - p_h_t||
    if isfield(to_log,'t') && isfield(to_log,'e_h') && ~isempty(to_log.t) && ~isempty(to_log.e_h)
        plot(to_log.t,  to_log.e_h,  'b-', 'LineWidth',1.2);
    end
    if isfield(sim_log,'t') && isfield(sim_log,'e_h') && ~isempty(sim_log.t) && ~isempty(sim_log.e_h)
        plot(sim_log.t, sim_log.e_h, 'r-', 'LineWidth',1.4);
    end

    % (b) hitch tracking error vs ref: ||p_h_r - p_h_ref||
    if isfield(ref,'t') && ~isempty(ref.t) ...
            && isfield(to_log,'t') && ~isempty(to_log.t) && isfield(to_log,'p_h_r') && ~isempty(to_log.p_h_r)
        p_h_ref_to = interp1(ref.t(:), p_h_ref, to_log.t(:), 'linear', 'extrap');
        e_track_to = sqrt(sum((to_log.p_h_r - p_h_ref_to).^2, 2));
        plot(to_log.t, e_track_to, 'b--', 'LineWidth',1.0);
    end
    if isfield(ref,'t') && ~isempty(ref.t) ...
            && isfield(sim_log,'t') && ~isempty(sim_log.t) && isfield(sim_log,'p_h_r') && ~isempty(sim_log.p_h_r)
        p_h_ref_sim = interp1(ref.t(:), p_h_ref, sim_log.t(:), 'linear', 'extrap');
        e_track_sim = sqrt(sum((sim_log.p_h_r - p_h_ref_sim).^2, 2));
        plot(sim_log.t, e_track_sim, 'r--', 'LineWidth',1.0);
    end

    legend({'to hitch overlap (constraint)', 'sim hitch overlap (constraint)', ...
            'to hitch tracking (vs ref)',     'sim hitch tracking (vs ref)'}, ...
            'Location','best');

    % ---------------- figure 4: qp hitch monitors (h/hdot/hddot) ----------------
    % expects (optional) fields:
    %   sim_log.qp_h [N x 2], sim_log.qp_hdot [N x 2], sim_log.qp_e_yddot [N x 2] (or sim_log.qp_yddot)
    %   to_log.qp_h  [N x 2], to_log.qp_hdot  [N x 2], to_log.qp_e_yddot  [N x 2] (or to_log.qp_yddot)

    has_sim_h    = isfield(sim_log,'qp_h')       && ~isempty(sim_log.qp_h);
    has_sim_hdot = isfield(sim_log,'qp_hdot')    && ~isempty(sim_log.qp_hdot);
    has_sim_ea   = isfield(sim_log,'qp_e_yddot') && ~isempty(sim_log.qp_e_yddot);
    has_sim_ya   = isfield(sim_log,'qp_yddot')   && ~isempty(sim_log.qp_yddot);

    has_to_h     = isfield(to_log,'qp_h')        && ~isempty(to_log.qp_h);
    has_to_hdot  = isfield(to_log,'qp_hdot')     && ~isempty(to_log.qp_hdot);
    has_to_ea    = isfield(to_log,'qp_e_yddot')  && ~isempty(to_log.qp_e_yddot);
    has_to_ya    = isfield(to_log,'qp_yddot')    && ~isempty(to_log.qp_yddot);

    if has_sim_h || has_sim_hdot || has_sim_ea || has_sim_ya || has_to_h || has_to_hdot || has_to_ea || has_to_ya
        figure('Color','w');
        tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

        % (1) ||h||
        nexttile; grid on; hold on;
        title('qp monitor: ||h|| (h = y - y_d)');
        xlabel('t (s)'); ylabel('||h||');
        if has_to_h && isfield(to_log,'t') && ~isempty(to_log.t)
            plot(to_log.t, sqrt(sum(to_log.qp_h.^2,2)), 'b-', 'LineWidth',1.1);
        end
        if has_sim_h && isfield(sim_log,'t') && ~isempty(sim_log.t)
            plot(sim_log.t, sqrt(sum(sim_log.qp_h.^2,2)), 'r-', 'LineWidth',1.3);
        end
        le = {};
        if has_to_h,  le{end+1}  = 'to ||h||';  end
        if has_sim_h, le{end+1}  = 'sim ||h||'; end
        if ~isempty(le), legend(le,'Location','best'); end

        % (2) ||hdot||
        nexttile; grid on; hold on;
        title('qp monitor: ||hdot|| (hdot = ydot - ydot_d)');
        xlabel('t (s)'); ylabel('||hdot||');
        if has_to_hdot && isfield(to_log,'t') && ~isempty(to_log.t)
            plot(to_log.t, sqrt(sum(to_log.qp_hdot.^2,2)), 'b-', 'LineWidth',1.1);
        end
        if has_sim_hdot && isfield(sim_log,'t') && ~isempty(sim_log.t)
            plot(sim_log.t, sqrt(sum(sim_log.qp_hdot.^2,2)), 'r-', 'LineWidth',1.3);
        end
        le = {};
        if has_to_hdot,  le{end+1}  = 'to ||hdot||';  end
        if has_sim_hdot, le{end+1}  = 'sim ||hdot||'; end
        if ~isempty(le), legend(le,'Location','best'); end

        % (3) ||e_yddot|| preferred, else ||yddot||
        nexttile; grid on; hold on;
        if (has_to_ea || has_sim_ea)
            title('qp monitor: ||e_{yddot}|| (e = yddot - yddot_{des})');
            ylabel('||e_{yddot}||');
            if has_to_ea && isfield(to_log,'t') && ~isempty(to_log.t)
                plot(to_log.t, sqrt(sum(to_log.qp_e_yddot.^2,2)), 'b-', 'LineWidth',1.1);
            end
            if has_sim_ea && isfield(sim_log,'t') && ~isempty(sim_log.t)
                plot(sim_log.t, sqrt(sum(sim_log.qp_e_yddot.^2,2)), 'r-', 'LineWidth',1.3);
            end
            le = {};
            if has_to_ea,  le{end+1}  = 'to ||e_{yddot}||';  end
            if has_sim_ea, le{end+1}  = 'sim ||e_{yddot}||'; end
        else
            title('qp monitor: ||yddot|| (actual)');
            ylabel('||yddot||');
            if has_to_ya && isfield(to_log,'t') && ~isempty(to_log.t)
                plot(to_log.t, sqrt(sum(to_log.qp_yddot.^2,2)), 'b-', 'LineWidth',1.1);
            end
            if has_sim_ya && isfield(sim_log,'t') && ~isempty(sim_log.t)
                plot(sim_log.t, sqrt(sum(sim_log.qp_yddot.^2,2)), 'r-', 'LineWidth',1.3);
            end
            le = {};
            if has_to_ya,  le{end+1}  = 'to ||yddot||';  end
            if has_sim_ya, le{end+1}  = 'sim ||yddot||'; end
        end
        xlabel('t (s)');
        if ~isempty(le), legend(le,'Location','best'); end
    end

    % ---------------- figure 5: input margin to bounds (sim) ----------------
    if isfield(sim_log,'u_margin') && ~isempty(sim_log.u_margin) && isfield(sim_log,'t') && ~isempty(sim_log.t)
        figure('Color','w'); grid on; hold on;
        title('input margin to bounds (sim)');
        xlabel('t (s)'); ylabel('min distance to bound');
        plot(sim_log.t, sim_log.u_margin(:,1), 'r-', 'LineWidth',1.2);
        plot(sim_log.t, sim_log.u_margin(:,2), 'b-', 'LineWidth',1.2);
        legend({'F margin','tau margin'},'Location','best');
    end
    fprintf('[sim] sat rate: F=%.1f%%, tau=%.1f%%\n', 100*mean(sim_log.u_sat(:,1)), 100*mean(sim_log.u_sat(:,2)));


    % ---------------- to cost breakdown (optional) ----------------
    if ~isempty(to_dbg) && isstruct(to_dbg) && isfield(to_dbg,'cost') && ~isempty(to_dbg.cost)
        J = to_dbg.cost;

        if isfield(to_dbg,'cost') && isfield(to_dbg.cost,'sum')
            cs = to_dbg.cost.sum;
            if isfield(cs,'pos') && isfield(cs,'u') && isfield(cs,'yaw') && isfield(cs,'h') && isfield(cs,'sm') && isfield(cs,'total_like')
                fprintf('cost sum: pos=%.3e, u=%.3e, yaw=%.3e, hitch=%.3e, sm=%.3e, total_like=%.3e\n', ...
                    cs.pos, cs.u, cs.yaw, cs.h, cs.sm, cs.total_like);
            end
        end

        if isfield(J,'J_pos') || isfield(J,'J_u') || isfield(J,'J_yaw') || isfield(J,'J_h')
            figure('Color','w'); hold on; grid on;
            if isfield(J,'J_pos'), plot(J.J_pos, 'LineWidth', 1); end
            if isfield(J,'J_u'),   plot(J.J_u,   'LineWidth', 1); end
            if isfield(J,'J_yaw'), plot(J.J_yaw, 'LineWidth', 1); end
            if isfield(J,'J_h'),   plot(J.J_h,   'LineWidth', 1); end
            xlabel('k'); ylabel('per-step cost');
            legend('J\_pos','J\_u','J\_yaw','J\_h');
            title('to cost breakdown (per step)');
        end

        if isfield(J,'J_vel') && ~isempty(J.J_vel)
            figure('Color','w'); grid on;
            plot(J.J_vel, 'LineWidth', 1);
            xlabel('k'); ylabel('J\_vel (monitor)');
            title('hitch velocity mismatch monitor');
        end

        if isfield(J,'J_pos') && ~isempty(J.J_pos)
            [~, kmax] = max(J.J_pos);
            if isfield(to_dbg,'dt') && ~isempty(to_dbg.dt)
                fprintf('worst J_pos at k=%d (t=%.3f s)\n', kmax, (kmax-1)*to_dbg.dt);
            else
                fprintf('worst J_pos at k=%d\n', kmax);
            end
        end
    end
    u = to_log.U;                      % N x 2
    F = u(:,1);  tau = u(:,2);
    
    figure; plot(abs(F)); hold on; yline(params.Fmax,'--'); grid on;
    title('|F| and bound'); xlabel('k');
    
    figure; plot(abs(tau)); hold on; yline(params.Taumax,'--'); grid on;
    title('|tau| and bound'); xlabel('k');
    
    satF  = find(abs(F)  > 0.99*params.Fmax);
    satTau= find(abs(tau)> 0.99*params.Taumax);
    fprintf('F near-sat steps: %d, tau near-sat steps: %d\n', numel(satF), numel(satTau));

    % ---------------- optional: analyze A_y conditioning ----------------
if isfield(sim_log,'A_y') && ~isempty(sim_log.A_y)

    Ay = sim_log.A_y;                 % expected: [N x 2 x 2]
    t  = sim_log.t(:);

    % handle the case Ay is [2 x 2 x N]
    if ndims(Ay) == 3 && size(Ay,1)==2 && size(Ay,2)==2 && size(Ay,3)==numel(t)
        Ay = permute(Ay, [3 1 2]);     % -> [N x 2 x 2]
    end

    N = size(Ay,1);
    if numel(t) ~= N
        error('A_y length mismatch: numel(t)=%d, size(A_y,1)=%d', numel(t), N);
    end

    condAy = nan(N,1);
    smax   = nan(N,1);
    smin   = nan(N,1);

    for k = 1:N
        A = squeeze(Ay(k,:,:));
        if any(~isfinite(A(:))), continue; end
        sv = svd(A);
        smax(k)   = sv(1);
        smin(k)   = sv(end);
        condAy(k) = sv(1) / max(sv(end), 1e-12);
    end

    fprintf('\n[sim] A_y stats:\n');
    fprintf('  cond(A_y): min=%.3g, median=%.3g, max=%.3g\n', ...
        min(condAy,[],'omitnan'), median(condAy,'omitnan'), max(condAy,[],'omitnan'));
    fprintf('  sigma_min(A_y): min=%.3g, median=%.3g\n', ...
        min(smin,[],'omitnan'), median(smin,'omitnan'));

    figure('Color','w'); grid on; hold on;
    title('A_y conditioning (sim)');
    xlabel('t (s)'); ylabel('cond(A_y)');
    plot(t, condAy, 'LineWidth', 1.2);

    figure('Color','w'); grid on; hold on;
    title('A_y singular values (sim)');
    xlabel('t (s)'); ylabel('singular values');
    plot(t, smax, 'LineWidth', 1.2);
    plot(t, smin, 'LineWidth', 1.2);
    legend({'sigma_max','sigma_min'}, 'Location', 'best');
end


    % ---------------- optional: rv if available ----------------
    if isfield(sim_log,'rv') && ~isempty(sim_log.rv) && isfield(sim_log,'t') && ~isempty(sim_log.t)
        figure('Color','w'); grid on; hold on;
        title('constraint velocity residual  ||Jv|| (sim)');
        xlabel('t (s)'); ylabel('rv');
        plot(sim_log.t, sim_log.rv, 'r-', 'LineWidth',1.3);
    end
end
