function [F_drive, tau_r] = fb_lin_controller(x, t, params, y_d, ydot_d, yddot_d)
% Full nonlinear feedback-linearizing controller for robot+trailer.
% Inputs:
%   x : 12x1 state [xr;yr;thetar; vxr;vyr;wr; xt;yt;thetat; vxt;vyt;wt]
%   t : time (unused here but kept for interface)
%   params: struct (m_r,I_r,m_t,I_t,d,Lt,...)
%   y_d, ydot_d, yddot_d: desired output trajectory (2x1 each)
% Outputs:
%   F_drive, tau_r : actuator commands
%
% Notes: This function recomputes M,J,J_dot etc each call (fine for sim).

% Unpack state
xr = x(1); yr = x(2); thetar = x(3);
vxr = x(4); vyr = x(5); wr = x(6);
xt = x(7); yt = x(8); thetat = x(9);
vxt = x(10); vyt = x(11); wt = x(12);

v = [vxr; vyr; wr; vxt; vyt; wt];

% Params
m_r = params.m_r; I_r = params.I_r;
m_t = params.m_t; I_t = params.I_t;
d   = params.d;  Lt = params.Lt;

B_q = [ cos(thetar), 0;
        sin(thetar), 0;
        0,           1;
        0,           0;
        0,           0;
        0,           0 ];

% Output Jacobian
J_y = [ 1, 0,  d*sin(thetar),  0, 0, 0;
        0, 1, -d*cos(thetar),  0, 0, 0 ];

% time derivative of J_y times qdot (Jdot_y * v)
% d/dt [ d*sin(thr) ] = d*cos(thr)*wr
Jdot_y_qdot = [ 0; 0 ];
Jdot_y_qdot(1) = d*cos(thetar)*wr * 1;
Jdot_y_qdot(2) = d*sin(thetar)*wr * 1;
Jdot_y_qdot = [ d*cos(thetar)*wr;
                 d*sin(thetar)*wr ];

r_rh_body = [-d; 0];           % robot hitch in robot body frame
r_th_body = [ Lt/2; 0];        % trailer hitch in trailer body frame

Rr = [cos(thetar) -sin(thetar); sin(thetar) cos(thetar)];
Rt = [cos(thetat) -sin(thetat); sin(thetat) cos(thetat)];

p_r = Rr * r_rh_body;   % world vector from robot COM to hitch
p_t = Rt * r_th_body;   % world vector from trailer COM to hitch

S_p_r = [-p_r(2); p_r(1)];
S_p_t = [-p_t(2); p_t(1)];

J_holo = [ eye(2), S_p_r, -eye(2), -S_p_t ];

% J_nonholo: l = [-sin(thetat); cos(thetat)]
l = [-sin(thetat); cos(thetat)];
J_nonholo = [0 0 0 l(1) l(2) 0];

% Full J
J = [ J_holo; J_nonholo ];   % 3x6

% S_p_r = [ d*sin(thr); -d*cos(thr) ]
Sd_p_r = [ d*cos(thetar)*wr; d*sin(thetar)*wr ]; 

% S_p_t = [-Lt/2*sin(tht); Lt/2*cos(tht)]
Sd_p_t = [ -(Lt/2)*cos(thetat)*wt; -(Lt/2)*sin(thetat)*wt ];

% Build dotJ_holo and dotJ_nonholo (each is 2x6 and 1x6 respectively)
dotJ_holo = [ zeros(2,2), Sd_p_r, zeros(2,2), -Sd_p_t ];
% dotJ_holo block layout: [0_2x2 , dS_p_r , 0_2x2 , -dS_p_t]

% dot of l = [-cos(tht)*wt; -sin(tht)*wt]
ldot = [-cos(thetat)*wt; -sin(thetat)*wt];
dotJ_nonholo = [ 0 0 0 ldot(1) ldot(2) 0 ];

% Full dotJ
dotJ = [ dotJ_holo; dotJ_nonholo ];  % 3x6

% Build M and H
M = diag([m_r m_r I_r m_t m_t I_t]);

H = zeros(6,1);
W = J * (M\J');      % 3x3
% numerical safeguard
invW = pinv(W);      % use pseudo-inverse for numerical robustness

% compute projection  
Minv = inv(M);   % safe since M diagonal
P = Minv - Minv * J' * (invW * (J * Minv));

% compute f and G = P*B_q
% f = P*H - M^{-1} J' * (invW * (dotJ * v))
f = P * H - Minv * J' * ( invW * (dotJ * v) );

G = P * B_q;   % 6x2

% compute Lf^2 y and Lg^2 y 
Lf2y = J_y * f + Jdot_y_qdot;   % 2x1
Lg2y = J_y * G;                 % 2x2

% desired linearized accel v (PD)
if nargin < 4 || isempty(y_d)
    y_d = [xr - d*cos(thetar); yr - d*sin(thetar)];  % hold current by default
end
if nargin < 5 || isempty(ydot_d)
    ydot_d = zeros(2,1);
end
if nargin < 6 || isempty(yddot_d)
    yddot_d = zeros(2,1);
end

% Gains (tune me)
Kp = diag([10, 10]);
Kd = diag([5, 5]);

% current output and derivative
y = [ xr - d*cos(thetar); yr - d*sin(thetar) ];
ydot = J_y * v;

v_des = yddot_d - Kd*(ydot - ydot_d) - Kp*(y - y_d);  % 2x1

% solve for u with damping regularization if needed 
eps = 1e-6;
if rcond(Lg2y) < 1e-6
    u = (Lg2y'*Lg2y + eps*eye(2)) \ (Lg2y' * (v_des - Lf2y));
else
    u = Lg2y \ (v_des - Lf2y);
end

% Saturation
if isfield(params,'Fmax'), Fmax = params.Fmax; else Fmax = 1e6; end
if isfield(params,'Taumax'), Taumax = params.Taumax; else Taumax = 1e6; end

F_drive = max(min(u(1), Fmax), -Fmax);
tau_r   = max(min(u(2), Taumax), -Taumax);
end
