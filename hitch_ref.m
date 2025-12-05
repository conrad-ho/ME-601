function [yd, ydot_d, yddot_ff, ref_state] = hitch_ref(t, y, ref_state)
% hitch_ref  generates reference trajectories for the hitch (y_d, ydot_d, yddot_ff)
%
% Inputs:
%   t          current time (scalar)
%   y          current hitch position [2x1]
%   ref_state  reference-trajectory state struct:
%              must include field:
%                .mode  = 'hold' | 'step_x' | 'local_circle' | 'L_turn'
%              optional fields depending on mode:
%                .y0       initial position
%                .dx, .Tf  (parameters for step_x)
%                .R, .w    (parameters for circle)
%                .v_ref, .T1, .R_turn (parameters for L_turn)
%
% Outputs:
%   yd, ydot_d, yddot_ff : desired position / velocity / feedforward accel
%   ref_state            : updated internal state (e.g., writes y0, initialization flag, etc.)

    % -------- initialize ref_state --------
    if nargin < 3 || isempty(ref_state)
        ref_state = struct();
        ref_state.mode = 'hold';
    end

    if ~isfield(ref_state,'initialized') || ~ref_state.initialized
        % starting point: only use input y on the first call
        if ~isfield(ref_state,'y0') || isempty(ref_state.y0)
            ref_state.y0 = y;      % afterwards always use y0, not y
        end
        % starting time: only use input t on the first call
        if ~isfield(ref_state,'t0') || isempty(ref_state.t0)
            ref_state.t0 = t;      % afterwards always use t_rel = t - t0
        end

        % default parameters for each mode
        switch ref_state.mode
            case 'step_x'
                if ~isfield(ref_state,'dx'),  ref_state.dx  = 0.2; end
                if ~isfield(ref_state,'Tf'),  ref_state.Tf  = 8.0; end

            case 'local_circle'
                if ~isfield(ref_state,'R'),   ref_state.R   = 1.0; end   % radius
                if ~isfield(ref_state,'w'),   ref_state.w   = 0.275; end % angular velocity
                % precompute circle geometry information
                R = ref_state.R;
                ref_state.phi0       = -pi/2;                 % starting tangent along +x
                ref_state.center     = ref_state.y0 + [0; R]; % center above y0 by R
                ref_state.angle_max  = 2*pi;                  % total swept angle (e.g. 90°)

            case 'L_turn'
                if ~isfield(ref_state,'v_ref'),   ref_state.v_ref   = 0.1; end
                if ~isfield(ref_state,'T1'),      ref_state.T1      = 8.0; end
                if ~isfield(ref_state,'R_turn'),  ref_state.R_turn  = 1.0; end

            case 'hold'
                % no additional params
        end

        ref_state.initialized = true;
    end

    % convenience aliases: from now on, only use y0 and t_rel
    y0    = ref_state.y0;
    t_rel = t - ref_state.t0;

    % =====================================================
    % generate trajectories based on mode
    % =====================================================
    switch ref_state.mode

        case 'hold'
            % stay at initial point
            yd       = y0;
            ydot_d   = [0; 0];
            yddot_ff = [0; 0];

        case 'step_x'
            % ------- S-curve that moves dx forward along +x direction -------
            dx = ref_state.dx;
            Tf = ref_state.Tf;

            if t_rel <= 0
                s = 0; s_dot = 0; s_ddot = 0;
            elseif t_rel >= Tf
                s = 1; s_dot = 0; s_ddot = 0;
            else
                tau = t_rel / Tf;         % in (0,1)
                s      = 3*tau^2 - 2*tau^3;
                s_dot  = (6*tau - 6*tau^2) / Tf;
                s_ddot = (6 - 12*tau) / Tf^2;
            end

            yd = y0 + [dx * s;
                       0];

            ydot_d = [dx * s_dot;
                      0];

            yddot_ff = [dx * s_ddot;
                        0];

        case 'local_circle'
            % smooth circular arc, starting at y0, initial tangent along +x,
            % circle center located above y0 by R
            R         = ref_state.R;
            w         = ref_state.w;          
            phi0      = ref_state.phi0;       
            center    = ref_state.center;     
            angle_max = ref_state.angle_max;  

            % total time for the circular arc
            T_arc = angle_max / abs(w);

            if t_rel <= 0
                % not started yet
                theta     = phi0;
                theta_dot = 0;
            elseif t_rel <= T_arc
                % moving along arc
                theta     = phi0 + w * t_rel;
                theta_dot = w;
            else
                % arc completed, hold endpoint
                theta     = phi0 + sign(w)*angle_max;
                theta_dot = 0;
            end

            % position
            yd = center + R*[cos(theta); sin(theta)];

            % velocity: y' = R*theta_dot*[-sinθ; cosθ]
            ydot_d = R * theta_dot * [-sin(theta); cos(theta)];

            % acceleration: for theta_ddot = 0, y'' = -R*theta_dot^2*[cosθ; sinθ]
            yddot_ff = -R * theta_dot^2 * [cos(theta); sin(theta)];

        case 'L_turn'
            % ------- straight line + left 90° turn -------
            v_ref  = ref_state.v_ref;
            T1     = ref_state.T1;
            R_turn = ref_state.R_turn;

            w      = v_ref / R_turn;
            T2     = (pi/2) / w;   % 90-degree turning duration

            line_end = y0 + [v_ref * T1; 0];
            center   = line_end + [0; -R_turn];  % left turn: center lies below

            if t_rel <= 0
                yd       = y0;
                ydot_d   = [0;0];
                yddot_ff = [0;0];

            elseif t_rel <= T1
                % segment 1: straight line along +x
                yd = y0 + [v_ref * t_rel;
                           0];

                ydot_d   = [v_ref; 0];
                yddot_ff = [0; 0];

            elseif t_rel <= T1 + T2
                % segment 2: quarter circle left turn
                tau   = t_rel - T1;
                theta = w * tau;

                yd = center + [ R_turn * sin(theta);
                                R_turn * (1 - cos(theta)) ];

                ydot_d = [  R_turn * w * cos(theta);
                            R_turn * w * sin(theta) ];

                yddot_ff = [ -R_turn * w^2 * sin(theta);
                              R_turn * w^2 * cos(theta) ];

            else
                % segment 3: hold at final point
                theta = w * T2;  % = pi/2
                yd = center + [ R_turn * sin(theta);
                                R_turn * (1 - cos(theta)) ];

                ydot_d   = [0; 0];
                yddot_ff = [0; 0];
            end

        otherwise
            error('Unknown ref_state.mode = %s', ref_state.mode);
    end
end
