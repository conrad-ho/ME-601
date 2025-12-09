function check_hitch_closure(x_to, params)
    d  = params.d;
    Lt = params.Lt;

    N = size(x_to,1);
    errs = zeros(N,1);

    for k = 1:N
        xk = x_to(k,:).';

        xr = xk(1); yr = xk(2); thetar = xk(3);
        xt = xk(7); yt = xk(8); thetat = xk(9);

        Rr = [cos(thetar) -sin(thetar);
              sin(thetar)  cos(thetar)];
        Rt = [cos(thetat) -sin(thetat);
              sin(thetat)  cos(thetat)];

        r_rh = [-d; 0];         % robot → hitch
        r_th = [Lt/2; 0];       % trailer COM → hitch

        p_h_r = [xr; yr] + Rr*r_rh;   % hitch (from robot)
        p_h_t = [xt; yt] - Rt*r_th;   % hitch (from trailer)

        e = p_h_r - p_h_t;
        errs(k) = norm(e);
    end

    fprintf('max hitch closure error = %.3e m\n', max(errs));
end
