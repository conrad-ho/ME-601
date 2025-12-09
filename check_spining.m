data = load('log.mat');
log  = data.log;

t_log      = log.t(:);
theta_log  = log.X(:,3);        % 仿真中的 robot yaw
theta_unw  = unwrap(theta_log); % 解除 2π wrap 看真实累计转角

figure;
subplot(2,1,1);
plot(t_log, theta_log, 'b-'); grid on;
ylabel('\theta_r (rad)');
title('Robot yaw in simulation (raw)');

subplot(2,1,2);
plot(t_log, theta_unw, 'r-'); grid on;
ylabel('unwrap(\theta_r) (rad)');
xlabel('time (s)');
title('Robot yaw in simulation (unwrapped)');
