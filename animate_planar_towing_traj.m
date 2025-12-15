function animate_planar_towing_traj(t, X, params, opts)
% animate_planar_towing_traj
% animate robot + trailer rectangles and their XY paths
%
% inputs:
%   t     : [N x 1] time
%   X     : [N x 12] state trajectory
%   params: struct with fields Wr,Lr,Wt,Lt
%   opts  : optional struct
%       .stride (default 5)    frame skip
%       .follow (default true) camera follow robot
%       .axis_pad (default [3 5 3 3]) [back front down up] meters around robot
%       .xlim (default [])     fixed axis if provided
%       .ylim (default [])     fixed axis if provided
%
% note: comments intentionally lower-case (per your preference)

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'stride'),   opts.stride = 5; end
    if ~isfield(opts,'follow'),   opts.follow = true; end
    if ~isfield(opts,'axis_pad'), opts.axis_pad = [3 5 3 3]; end
    if ~isfield(opts,'xlim'),     opts.xlim = []; end
    if ~isfield(opts,'ylim'),     opts.ylim = []; end

    % figure setup
    figure('Color','w'); hold on; grid on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)');
    title('Planar Robot–Trailer Dynamics (Rigid Hitch)');

    if ~isempty(opts.xlim), xlim(opts.xlim); else, xlim([-2 12]); end
    if ~isempty(opts.ylim), ylim(opts.ylim); else, ylim([-4 4]); end

    path_robot   = plot(NaN,NaN,'k-','LineWidth',1.2,'DisplayName','robot');
    path_trailer = plot(NaN,NaN,'k--','LineWidth',1.2,'DisplayName','trailer');

    robot_patch   = patch(NaN,NaN,'r','FaceAlpha',0.4,'EdgeColor','none');
    trailer_patch = patch(NaN,NaN,'b','FaceAlpha',0.4,'EdgeColor','none');

    legend('Location','best');

    idx_vec = 1:opts.stride:length(t);
    Nsteps  = numel(idx_vec);

    traj_x_r = NaN(1, Nsteps); traj_y_r = NaN(1, Nsteps);
    traj_x_t = NaN(1, Nsteps); traj_y_t = NaN(1, Nsteps);

    kplot = 0;
    for ii = 1:Nsteps
        i = idx_vec(ii);
        kplot = kplot + 1;

        xr = X(i,1); yr = X(i,2); thetar = X(i,3);
        xt = X(i,7); yt = X(i,8); thetat = X(i,9);

        traj_x_r(kplot) = xr; traj_y_r(kplot) = yr;
        traj_x_t(kplot) = xt; traj_y_t(kplot) = yt;

        set(path_robot,  'XData', traj_x_r(1:kplot), 'YData', traj_y_r(1:kplot));
        set(path_trailer,'XData', traj_x_t(1:kplot), 'YData', traj_y_t(1:kplot));

        set(robot_patch,'XData',rect_x(xr,params.Wr,params.Lr,thetar), ...
                       'YData',rect_y(yr,params.Wr,params.Lr,thetar));
        set(trailer_patch,'XData',rect_x(xt,params.Wt,params.Lt,thetat), ...
                         'YData',rect_y(yt,params.Wt,params.Lt,thetat));

        if opts.follow
            pad = opts.axis_pad;
            axis([xr-pad(1) xr+pad(2) yr-pad(3) yr+pad(4)]);
        end

        drawnow;
    end
end

% local helpers (kept identical to your original)
function Xpoly = rect_x(xc,W,L,theta)
    R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
    pts=0.5*[-L -W; L -W; L W; -L W]';
    Xpoly = R(1,:)*pts + xc;
end

function Ypoly = rect_y(yc,W,L,theta)
    R=[cos(theta) -sin(theta); sin(theta) cos(theta)];
    pts=0.5*[-L -W; L -W; L W; -L W]';
    Ypoly = R(2,:)*pts + yc;
end
