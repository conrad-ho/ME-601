function check_to_states_vs_dynamics(x_to, params)
    N = size(x_to, 2);
    bad_idx = [];

    for k = 1:N
        xk = x_to(:,k);
        try
            xdot_k = towing_dynamics_full(xk, [0;0], params);
        catch ME
            fprintf(2, 'towing_dynamics_full crashed at TO node k=%d\n', k);
            fprintf(2, '  message: %s\n', ME.message);
            bad_idx(end+1) = k; %#ok<AGROW>
            continue;
        end
    end

    if isempty(bad_idx)
        fprintf('All TO states are acceptable for towing_dynamics_full (u=0).\n');
    else
        fprintf(2, 'TO states cause dynamics NaN at indices: ');
        fprintf(2, '%d ', bad_idx);
        fprintf(2, '\n');
    end
end
