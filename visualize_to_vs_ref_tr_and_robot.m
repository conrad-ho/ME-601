function visualize_to_vs_ref_tr_and_robot(x_to, ref_to, params, dt_to)
%VISUALIZE_TO_VS_REF
%   Visualize TO results vs reference paths (robot + trailer)
%
% x_to   : [nx x (N+1)] or (N+1 x nx) state trajectory from TO
%          state = [xr, yr, thetar, vxr, vyr, wr,  xt, yt, thetat, vxt, vyt, wt]
%
% ref_to : struct containing:
%           .t    [1 x (N+1)]      time vector
%           .p_r  [2 x (N+1)]      reference robot COM path
%           .p_tr [2 x (N+1)]      reference trailer COM path
%
% params : system parameter struct (not used here but kept for interface)
% dt_to  : nominal time step, used if ref_to.t is unavailable

    % ---------- Keep your original transpose behavior ----------
    x_to = x_to.';
    [nx, Ncol] = size(x_to); %#ok<ASGLU>

    % ---------- Time vector ----------
    if isfield(ref_to, 't') && numel(ref_to.t) == Ncol
        t = ref_to.t(:)';   % 1 x (N+1)
    else
        t = (0:Ncol-1) * dt_to;
    end

    % ---------- Extract robot / trailer COM positions ----------
    p_r_to  = zeros(2, Ncol);
    p_tr_to = zeros(2, Ncol);

    for k = 1:Ncol
        xk = x_to(:,k);

        % x = [xr, yr, thetar, vxr, vyr, wr,  xt, yt, thetat, vxt, vyt, wt]
        p_r_to(:,k)  = [xk(1); xk(2)];
        p_tr_to(:,k) = [xk(7); xk(8)];
    end

    % ---------- Compute tracking errors ----------
    e_tr = [];
    e_r  = [];

    if isfield(ref_to, 'p_tr')
        e_tr = p_tr_to - ref_to.p_tr;
    end
    if isfield(ref_to, 'p_r')
        e_r = p_r_to - ref_to.p_r;
    end

    % =====================================================================
    % Figure 1: XY trajectory comparison
    % =====================================================================
    figure;
    tiledlayout(1,2);

    % --- (a) XY trajectories ---
    nexttile;
    hold on; grid on; axis equal;

    % Reference paths
    if isfield(ref_to, 'p_tr')
        plot(ref_to.p_tr(1,:), ref_to.p_tr(2,:), 'k--', 'LineWidth', 1.5);
    end
    if isfield(ref_to, 'p_r')
        plot(ref_to.p_r(1,:),  ref_to.p_r(2,:),  'k:',  'LineWidth', 1.5);
    end

    % TO results
    plot(p_tr_to(1,:), p_tr_to(2,:), 'b-', 'LineWidth', 1.5);
    plot(p_r_to(1,:),  p_r_to(2,:),  'r-', 'LineWidth', 1.5);

    % Legend
    legend_entries = {};
    if isfield(ref_to, 'p_tr'), legend_entries{end+1} = 'ref trailer'; end
    if isfield(ref_to, 'p_r'),  legend_entries{end+1} = 'ref robot';   end
    legend_entries{end+1} = 'TO trailer';
    legend_entries{end+1} = 'TO robot';
    legend(legend_entries, 'Location', 'best');

    xlabel('x (m)');
    ylabel('y (m)');
    title('XY trajectories: reference vs TO');

    % Mark start points
    plot(p_tr_to(1,1), p_tr_to(2,1), 'bo', 'MarkerFaceColor','b');
    plot(p_r_to(1,1),  p_r_to(2,1),  'ro', 'MarkerFaceColor','r');

    % =====================================================================
    % Figure 2: tracking error vs time
    % =====================================================================
    nexttile;
    hold on; grid on;

    legend_list = {};

    if ~isempty(e_tr)
        e_tr_norm = vecnorm(e_tr);
        plot(t, e_tr_norm, 'b-', 'LineWidth', 1.5);
        legend_list{end+1} = '|e_{trailer}|';
    end

    if ~isempty(e_r)
        e_r_norm = vecnorm(e_r);
        plot(t, e_r_norm, 'r--', 'LineWidth', 1.5);
        legend_list{end+1} = '|e_{robot}|';
    end

    xlabel('time (s)');
    ylabel('position error (m)');

    if ~isempty(legend_list)
        legend(legend_list, 'Location', 'best');
    end
    title('Tracking error vs time');

    % =====================================================================
    % Figure 3: phi = thetat - thetar
    % =====================================================================
    phi = x_to(9,:) - x_to(3,:);   % trailer yaw - robot yaw
    figure;
    plot(t, phi, 'LineWidth', 1.5);
    grid on;
    xlabel('time (s)');
    ylabel('\phi (rad)');
    title('\phi (trailer yaw - robot yaw)');

end
