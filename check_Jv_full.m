function check_Jv_full(x_to, params, dt)

    [Nplus1, nx] = size(x_to);
    N = Nplus1 - 1;

    e_Jv_full = zeros(Nplus1, 1);
    e_Jv_xy   = zeros(Nplus1, 1);
    e_Jv_yaw  = zeros(Nplus1, 1);

    % For each timestep:
    for k = 1:Nplus1
        xk = x_to(k,:)';

        % Extract velocities
        v = [
            xk(4);  % v_xr
            xk(5);  % v_yr
            xk(6);  % w_r
            xk(10); % v_xt
            xk(11); % v_yt
            xk(12)  % w_t
        ];

        % --- POSITION CONSTRAINTS (as before) ---
        dyn = towing_dynamics_mats(xk, params);
        J_xy = dyn.J;  % This is 2x6

        % --- YAW CONSISTENCY CONSTRAINT ---
        % θ_t - θ_r = 0
        J_yaw = [0 0 -1  0 0 1];

        % --- FULL J ---
        J_full = [J_xy; J_yaw];

        % compute violation
        e_Jv_full(k) = norm(J_full * v);
        e_Jv_xy(k)   = norm(J_xy   * v);
        e_Jv_yaw(k)  = abs(J_yaw * v);
    end

    fprintf("=== FULL velocity constraint check ===\n");
    fprintf("max ||J_full v|| = %.3e\n", max(e_Jv_full));
    fprintf("RMS ||J_full v|| = %.3e\n", rms(e_Jv_full));

    fprintf("\n=== XY constraint (your original) ===\n");
    fprintf("max ||J_xy v|| = %.3e\n", max(e_Jv_xy));

    fprintf("\n=== yaw constraint ===\n");
    fprintf("max |J_yaw v| = %.3e\n", max(e_Jv_yaw));

    % optional plot
    t = (0:N)*dt;
    figure; plot(t, e_Jv_full); title('FULL ||J v||');
end
