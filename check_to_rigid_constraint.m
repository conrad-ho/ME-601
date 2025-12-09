function check_to_rigid_constraint(x_to, params)
%CHECK_TO_RIGID_CONSTRAINT
%   Quick diagnostic for TO trajectory x_to:
%   - checks how well the rigid hitch position constraint is satisfied.
%
% Usage:
%   1) Full dynamics state (N+1 x 12):
%        x = [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
%      Hitch positions:
%        p_hr = [xr - d*cos(thetar);
%                yr - d*sin(thetar)];
%        p_ht = [xt + Lt*cos(thetat);
%                yt + Lt*sin(thetat)];   % (change +Lt to -Lt if your geometry differs)
%
%   2) Reduced kinematic state (N+1 x 4):
%        x = [xr, yr, thetar, theta_t]
%      We reconstruct trailer COM from hitch:
%        p_h  = p_hr
%        p_tr = p_h - [Lt*cos(theta_t); Lt*sin(theta_t)];
%      and then reconstruct hitch from trailer side:
%        p_ht = p_tr + [Lt*cos(theta_t); Lt*sin(theta_t)];
%
% params : struct with at least fields
%          .d  : distance from robot COM to hitch (along -x_r)
%          .Lt : distance from trailer COM to hitch (along +/-x_t)

    if ~isfield(params,'d') || ~isfield(params,'Lt')
        error('params must contain fields d and Lt.');
    end

    d  = params.d;
    Lt = params.Lt;

    [Np, nx] = size(x_to);

    if nx ~= 12 && nx ~= 4
        error('x_to must be (N+1) x 12 or (N+1) x 4.');
    end

    p_hr = zeros(2,Np);  % hitch from robot side
    p_ht = zeros(2,Np);  % hitch from trailer side
    e_h  = zeros(2,Np);  % position error

    for k = 1:Np
        xk = x_to(k,:).';

        if nx == 12
            % -------- full-state version --------
            % x = [xr, yr, thetar, vxr, vyr, wr, xt, yt, thetat, vxt, vyt, wt]
            xr     = xk(1);
            yr     = xk(2);
            thetar = xk(3);
            xt     = xk(7);
            yt     = xk(8);
            thetat = xk(9);

            % hitch position from robot
            p_hr(:,k) = [xr - d*cos(thetar);
                         yr - d*sin(thetar)];

            % hitch position from trailer (change +Lt to -Lt if your geometry differs)
            p_ht(:,k) = [xt + Lt*cos(thetat);
                         yt + Lt*sin(thetat)];

        elseif nx == 4
            % -------- reduced kinematic version --------
            % x = [xr, yr, thetar, theta_t]
            xr      = xk(1);
            yr      = xk(2);
            thetar  = xk(3);
            theta_t = xk(4);

            % hitch from robot COM
            p_hr(:,k) = [xr - d*cos(thetar);
                         yr - d*sin(thetar)];

            % trailer COM reconstructed from hitch (trailer on the left of hitch)
            p_tr = p_hr(:,k) - [Lt*cos(theta_t);
                                Lt*sin(theta_t)];

            % hitch from trailer side (should equal p_hr if geometry is consistent)
            % NOTE: if your trailer geometry uses -Lt instead of +Lt in full model,
            %       flip the signs here consistently.
            p_ht(:,k) = p_tr + [Lt*cos(theta_t);
                                Lt*sin(theta_t)];
        end

        e_h(:,k) = p_hr(:,k) - p_ht(:,k);
    end

    e_norm = vecnorm(e_h,2,1);

    fprintf('--- Rigid hitch constraint diagnostics (position level) ---\n');
    fprintf('  state dim nx        = %d\n', nx);
    fprintf('  max ||p_hr - p_ht|| = %.3e m\n', max(e_norm));
    fprintf('  RMS ||p_hr - p_ht|| = %.3e m\n', sqrt(mean(e_norm.^2)));
    fprintf('  start error         = [%.3e, %.3e] m\n', e_h(1,1), e_h(2,1));
    fprintf('  end   error         = [%.3e, %.3e] m\n', e_h(1,end), e_h(2,end));
end
