function ref_to = lift_kin_to_full_ref(x_kin, params, dt)
% lift kinematic TO trajectory [xr;yr;thetar;theta_t]
% into a ref_to struct compatible with your existing code:
%   ref_to.t    : 1 x (N+1)
%   ref_to.p_r  : 2 x (N+1)  robot COM path
%   ref_to.p_tr : 2 x (N+1)  trailer COM path

    d  = params.d;
    Lt = params.Lt;

    Np = size(x_kin,1);    % N+1
    N  = Np - 1;

    t  = (0:N) * dt;

    p_r  = zeros(2,Np);
    p_tr = zeros(2,Np);

    for k = 1:Np
        xr      = x_kin(k,1);
        yr      = x_kin(k,2);
        thetar  = x_kin(k,3);
        theta_t = x_kin(k,4);

        % robot COM
        p_r(:,k) = [xr; yr];

        % hitch from robot com
        p_h = [xr - d*cos(thetar);
               yr - d*sin(thetar)];

        % trailer com
        p_tr(:,k) = p_h - [Lt*cos(theta_t);
                           Lt*sin(theta_t)];
    end

    ref_to = struct;
    ref_to.t    = t;
    ref_to.p_r  = p_r;
    ref_to.p_tr = p_tr;
end
