function test_calculate_cycle()
% TEST_CALCULATE_CYCLE 验证超临界CO₂布雷顿循环计算程序的正确性
%
% 本程序包含以下验证内容：
% 1. 基本功能测试 - 使用典型参数调用函数
% 2. 物理合理性检查 - 验证能量平衡和物理参数的合理性
% 3. 参数敏感性测试 - 分析关键参数变化对系统性能的影响
% 4. T-s图绘制 - 可视化循环过程
% 5. 收敛性测试 - 验证迭代计算的稳定性

% 设置典型工况参数
params = struct();
params.P_high = 16;          % MPa (修改为小于熔化压力)
params.P_low = 7.5;          % MPa
params.T_high = 873.15;      % K (600°C)
params.T_low = 305.15;       % K (32°C)
params.P_reheat = 10;        % MPa (修改为合适的再热压力)
params.P_intercool = 10;     % MPa (修改为合适的中间冷却压力)
params.split_ratio = 0.3;    % 分流比例
params.HT_recuperator_dT = 10; % K
params.LT_recuperator_dT = 5;  % K
params.eta_turbine_HP = 0.90;  % 高压透平效率
params.eta_turbine_LP = 0.90;  % 低压透平效率
params.eta_compressor_main = 0.89;  % 主压缩机效率
params.eta_compressor_recomp = 0.89; % 副压缩机效率
params.eta_recuperator_HT = 0.95;   % 高温回热器效率
params.eta_recuperator_LT = 0.95;   % 低温回热器效率
params.eta_heater = 0.95;           % 加热器效率
params.m_dot = 100;                % kg/s (降低质量流量以提高数值稳定性)

fprintf('开始验证超临界CO₂布雷顿循环计算程序...\n\n');

% 1. 基本功能测试
fprintf('1. 基本功能测试\n');
fprintf('使用典型参数计算循环性能...\n');
[states, performance] = calculate_cycle(params);

% 输出主要性能指标
fprintf('循环热效率: %.2f%%\n', performance.eta_thermal * 100);
fprintf('净功率: %.2f MW\n', performance.W_net / 1000);
fprintf('总透平功率: %.2f MW\n', performance.W_turbine / 1000);
fprintf('总压缩功率: %.2f MW\n', performance.W_compressor / 1000);
fprintf('总热输入: %.2f MW\n', performance.Q_in / 1000);
fprintf('\n');

% 【改进】打印低温回热器出口温度和热负荷
fprintf('【改进】T5=%.1fK, T6=%.1fK, Q_recup_LT=%.1f kW\n', ...
        states(5).T, states(6).T, performance.Q_recup_LT);

% 2. 物理合理性检查
fprintf('2. 物理合理性检查\n');

% 2.1 能量平衡检查
energy_balance = abs(performance.Q_in - (performance.W_net + performance.Q_cooler + performance.Q_intercooler));
fprintf('能量平衡误差: %.2f kW (应接近0)\n', energy_balance);

% 2.2 温度检查
temp_valid = true;
for i = 1:17
    if states(i).T < params.T_low || states(i).T > params.T_high
        temp_valid = false;
        fprintf('警告: 状态点%d温度超出范围: %.2f K\n', i, states(i).T);
    end
end
if temp_valid
    fprintf('所有状态点温度在合理范围内\n');
end

% 2.3 压力检查
pressure_valid = true;
for i = 1:17
    if states(i).P < params.P_low || states(i).P > params.P_high
        pressure_valid = false;
        fprintf('警告: 状态点%d压力超出范围: %.2f MPa\n', i, states(i).P);
    end
end
if pressure_valid
    fprintf('所有状态点压力在合理范围内\n');
end

% 2.4 回热器有效性检查
HT_effectiveness = (states(16).T - states(15).T) / (states(4).T - states(15).T);
LT_effectiveness = (states(15).T - states(14).T) / (states(5).T - states(14).T);
fprintf('高温回热器有效性: %.2f (应小于1)\n', HT_effectiveness);
fprintf('低温回热器有效性: %.2f (应小于1)\n', LT_effectiveness);
fprintf('\n');

fprintf('【调试】Q_recup_HT = %.1f kW, Q_recup_LT = %.1f kW\n', performance.Q_recup_HT, performance.Q_recup_LT);
fprintf('【调试】T5=%.1f, T14=%.1f, T15=%.1f, T16=%.1f\n', states(5).T, states(14).T, states(15).T, states(16).T);

% 3. 参数敏感性测试
fprintf('3. 参数敏感性测试\n');

% 3.1 最高温度敏感性
fprintf('最高温度敏感性分析:\n');
temps = [823.15, 848.15, 873.15, 898.15, 923.15];  % K
etas = zeros(size(temps));
for i = 1:length(temps)
    params_temp = params;
    params_temp.T_high = temps(i);
    [~, perf] = calculate_cycle(params_temp);
    etas(i) = perf.eta_thermal;
end
fprintf('温度(°C)\t热效率(%%)\n');
for i = 1:length(temps)
    fprintf('%.1f\t\t%.2f\n', temps(i)-273.15, etas(i)*100);
end
fprintf('\n');

% 4. T-s图绘制
fprintf('4. 绘制T-s图\n');
figure('Name', '超临界CO₂布雷顿循环T-s图');
hold on;

% 提取所有状态点的T和s值
T_points = zeros(1, 17);
s_points = zeros(1, 17);
for i = 1:17
    T_points(i) = states(i).T;
    s_points(i) = states(i).s;
end

% 绘制状态点和连接线
plot(s_points, T_points, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 15);

% 添加状态点标签
for i = 1:17
    text(s_points(i), T_points(i), sprintf(' %d', i), 'FontSize', 10);
end

% 设置图形属性
xlabel('比熵 s (kJ/kg·K)');
ylabel('温度 T (K)');
title('超临界CO₂布雷顿循环T-s图');
grid on;
hold off;

% 5. 收敛性测试
fprintf('5. 收敛性测试\n');
fprintf('迭代次数: %d\n', performance.iter_count);
if performance.iter_count < 100
    fprintf('迭代计算正常收敛\n');
else
    fprintf('警告: 迭代计算可能存在收敛问题\n');
end

fprintf('\n验证完成。\n');
end 