% 参数初始化
P = [0.5, 0.5; 0.5, 0.5]; % 状态转移矩阵
Q = 0.01; % 状态噪声
R = 0.1; % 观测噪声
A_modes = [1, -1]; % 模型模式: 1 对应状态0, -1 对应状态1

n_steps = 100; % 仿真步数
state_true = zeros(1, n_steps); % 真实状态
state_est = zeros(1, n_steps); % 估计状态
z = zeros(1, n_steps); % 观测值

% 初始状态
current_state = 0; 

% 模拟系统跳变和观测
for k = 1:n_steps
    if rand < P(current_state+1, 2) % 状态跳变判决
        current_state = 1 - current_state;
    end
    state_true(k) = current_state;
    z(k) = A_modes(current_state+1) + sqrt(R) * randn(); % 观测值
end

% IMM 估计
% 初始化概率
mu = [0.5, 0.5]; 
x_est = [1; -1]; % 模型下的状态估计
P_est = [1, 1]; 

for k = 1:n_steps
    % 模型1估计
    x_pred_1 = A_modes(1) * x_est(1);
    P_pred_1 = P_est(1) + Q;
    K1 = P_pred_1 / (P_pred_1 + R);
    x_est_1 = x_pred_1 + K1 * (z(k) - x_pred_1);
    P_est_1 = (1 - K1) * P_pred_1;
    
    % 模型2估计
    x_pred_2 = A_modes(2) * x_est(2);
    P_pred_2 = P_est(2) + Q;
    K2 = P_pred_2 / (P_pred_2 + R);
    x_est_2 = x_pred_2 + K2 * (z(k) - x_pred_2);
    P_est_2 = (1 - K2) * P_pred_2;
    
    % 模型概率更新
    L1 = exp(-0.5 * (z(k) - x_pred_1)^2 / R);
    L2 = exp(-0.5 * (z(k) - x_pred_2)^2 / R);
    mu(1) = mu(1) * L1;
    mu(2) = mu(2) * L2;
    mu = mu / sum(mu); % 归一化
    
    % 综合估计
    state_est(k) = mu(1) * x_est_1 + mu(2) * x_est_2;
    
    % 更新
    x_est = [x_est_1; x_est_2];
    P_est = [P_est_1, P_est_2];
end

% 绘图
figure;
subplot(2,1,1);
stairs(1:n_steps, state_true, 'r', 'LineWidth', 1.5); hold on;
stairs(1:n_steps, state_est < 0, 'b--', 'LineWidth', 1.5);
stairs(1:n_steps, z < 0, 'c:', 'LineWidth', 1.5);
legend('真实状态', '估计状态');
xlabel('时间步'); ylabel('状态');
title('跳变参数追踪实现二元判决');

subplot(2,1,2);
plot(1:n_steps, z, 'k');
xlabel('时间步'); ylabel('观测值');
title('观测值');