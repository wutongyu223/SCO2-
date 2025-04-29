%% SCO2_ParameterAnalysis.m - 超临界CO2布雷顿循环参数分析程序
% 本程序用于分析不同参数对超临界CO2布雷顿循环性能的影响
% 基于calautatlion_renum.m开发，但专注于参数扫描和图表生成

%% 设置默认参数（基准工况）
% 透平参数
default_n_hp_turb = 0.93; % 高压透平效率
default_n_lp_turb = 0.93; % 低压透平效率
default_P1 = 21;          % 高压透平入口压力(MPa)
default_P2 = 15;          % 高压透平出口压力(MPa)
default_P4 = 7.5729;      % 低压透平出口压力(MPa)
default_T1 = 873;         % 高压透平入口温度(K)

% 压缩机参数
default_N_mc_a = 0.89;    % 主压缩机a效率
default_N_mc_b = 0.89;    % 主压缩机b效率
default_N_rc = 0.85;      % 副压缩机效率
default_P8 = 12;          % 主压缩机a出口压力(MPa)
default_T7 = 305;         % 冷却器出口温度/主压缩机入口温度(K)

% 回热器参数
default_e1 = 0.86;        % 低温回热器效率
default_e2 = 0.86;        % 高温回热器效率

% 分流比
default_alpha = 0.3333;   % 默认分流比

% 其他参数
entropy_gen_factor = 0.03; % 合流过程熵增系数

%% 图3.1 分流比对循环热效率的影响
disp('正在计算：分流比对循环热效率的影响...');
alpha_range = 0.1:0.02:0.5;
n_alpha = length(alpha_range);
eta_alpha = zeros(1, n_alpha);
w_net_alpha = zeros(1, n_alpha);

for i = 1:n_alpha
    [eta_alpha(i), w_net_alpha(i), ~, ~] = calculate_cycle_performance(...
        default_P1, default_P2, default_P4, default_P8, alpha_range(i), ...
        default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
        default_e1, default_e2, default_T1, default_T7);
end

figure(1);
yyaxis left;
plot(alpha_range, eta_alpha*100, 'b-o', 'LineWidth', 1.5);
xlabel('分流比');
ylabel('循环热效率 (%)');
grid on;

yyaxis right;
plot(alpha_range, w_net_alpha/1000, 'r--*', 'LineWidth', 1.5);
ylabel('净输出功率 (kJ/kg)');
title('图3.1 分流比对循环热效率和净功率的影响');
legend('循环热效率', '净输出功率', 'Location', 'best');
savefig('图3.1_分流比影响.fig');
saveas(gcf, '图3.1_分流比影响.png');

%% 图3.2 透平入口温度对循环效率的影响
disp('正在计算：透平入口温度对循环效率的影响...');
T1_range = 773:20:973;
n_T1 = length(T1_range);
eta_T1 = zeros(1, n_T1);
w_net_T1 = zeros(1, n_T1);

for i = 1:n_T1
    [eta_T1(i), w_net_T1(i), ~, ~] = calculate_cycle_performance(...
        default_P1, default_P2, default_P4, default_P8, default_alpha, ...
        default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
        default_e1, default_e2, T1_range(i), default_T7);
end

figure(2);
yyaxis left;
plot(T1_range, eta_T1*100, 'b-o', 'LineWidth', 1.5);
xlabel('透平入口温度 (K)');
ylabel('循环热效率 (%)');
grid on;

yyaxis right;
plot(T1_range, w_net_T1/1000, 'r--*', 'LineWidth', 1.5);
ylabel('净输出功率 (kJ/kg)');
title('图3.2 透平入口温度对循环效率和净功率的影响');
legend('循环热效率', '净输出功率', 'Location', 'best');
savefig('图3.2_透平入口温度影响.fig');
saveas(gcf, '图3.2_透平入口温度影响.png');

%% 图3.3 透平入口压力对循环效率的影响
disp('正在计算：透平入口压力对循环效率的影响...');
P1_range = 18:0.5:25;
n_P1 = length(P1_range);
eta_P1 = zeros(1, n_P1);
w_net_P1 = zeros(1, n_P1);

for i = 1:n_P1
    % 保持压力比例不变
    P2_ratio = default_P2 / default_P1;
    P2 = P1_range(i) * P2_ratio;
    
    [eta_P1(i), w_net_P1(i), ~, ~] = calculate_cycle_performance(...
        P1_range(i), P2, default_P4, default_P8, default_alpha, ...
        default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
        default_e1, default_e2, default_T1, default_T7);
end

figure(3);
yyaxis left;
plot(P1_range, eta_P1*100, 'b-o', 'LineWidth', 1.5);
xlabel('透平入口压力 (MPa)');
ylabel('循环热效率 (%)');
grid on;

yyaxis right;
plot(P1_range, w_net_P1/1000, 'r--*', 'LineWidth', 1.5);
ylabel('净输出功率 (kJ/kg)');
title('图3.3 透平入口压力对循环效率和净功率的影响');
legend('循环热效率', '净输出功率', 'Location', 'best');
savefig('图3.3_透平入口压力影响.fig');
saveas(gcf, '图3.3_透平入口压力影响.png');

%% 图3.4 主压缩机入口温度对循环效率的影响
disp('正在计算：主压缩机入口温度对循环效率的影响...');
T7_range = 300:1:315;
n_T7 = length(T7_range);
eta_T7 = zeros(1, n_T7);
w_net_T7 = zeros(1, n_T7);

for i = 1:n_T7
    [eta_T7(i), w_net_T7(i), ~, ~] = calculate_cycle_performance(...
        default_P1, default_P2, default_P4, default_P8, default_alpha, ...
        default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
        default_e1, default_e2, default_T1, T7_range(i));
end

figure(4);
yyaxis left;
plot(T7_range, eta_T7*100, 'b-o', 'LineWidth', 1.5);
xlabel('主压缩机入口温度 (K)');
ylabel('循环热效率 (%)');
grid on;

yyaxis right;
plot(T7_range, w_net_T7/1000, 'r--*', 'LineWidth', 1.5);
ylabel('净输出功率 (kJ/kg)');
title('图3.4 主压缩机入口温度对循环效率和净功率的影响');
legend('循环热效率', '净输出功率', 'Location', 'best');
savefig('图3.4_主压缩机入口温度影响.fig');
saveas(gcf, '图3.4_主压缩机入口温度影响.png');

%% 图3.5 主压缩机出口压力对循环效率的影响
disp('正在计算：主压缩机出口压力对循环效率的影响...');
P8_range = 9:0.4:15;
n_P8 = length(P8_range);
eta_P8 = zeros(1, n_P8);
w_net_P8 = zeros(1, n_P8);

for i = 1:n_P8
    [eta_P8(i), w_net_P8(i), ~, ~] = calculate_cycle_performance(...
        default_P1, default_P2, default_P4, P8_range(i), default_alpha, ...
        default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
        default_e1, default_e2, default_T1, default_T7);
end

figure(5);
yyaxis left;
plot(P8_range, eta_P8*100, 'b-o', 'LineWidth', 1.5);
xlabel('主压缩机出口压力 (MPa)');
ylabel('循环热效率 (%)');
grid on;

yyaxis right;
plot(P8_range, w_net_P8/1000, 'r--*', 'LineWidth', 1.5);
ylabel('净输出功率 (kJ/kg)');
title('图3.5 主压缩机出口压力对循环效率和净功率的影响');
legend('循环热效率', '净输出功率', 'Location', 'best');
savefig('图3.5_主压缩机出口压力影响.fig');
saveas(gcf, '图3.5_主压缩机出口压力影响.png');

%% 图3.6 系统各设备㶲损分析
disp('正在计算：系统各设备㶲损分析...');
[~, ~, exergy_destruction, exergy_efficiency] = calculate_cycle_performance(...
    default_P1, default_P2, default_P4, default_P8, default_alpha, ...
    default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
    default_e1, default_e2, default_T1, default_T7);

% 设备名称
devices = {'高压透平', '低压透平', '主压缩机a', '主压缩机b', '副压缩机', ...
           '高温回热器', '低温回热器', '冷却器', '合流点', '加热器'};

% 㶲损比例
exergy_destruction_percentage = exergy_destruction / sum(exergy_destruction) * 100;

% 绘制饼图
figure(6);
pie(exergy_destruction_percentage);
legend(devices, 'Location', 'eastoutside');
title('图3.6 系统各设备㶲损分布');
savefig('图3.6_系统各设备㶲损.fig');
saveas(gcf, '图3.6_系统各设备㶲损.png');

% 绘制柱状图
figure(7);
bar([exergy_destruction; exergy_efficiency]', 'grouped');
set(gca, 'XTickLabel', devices);
xtickangle(45);
grid on;
legend('㶲损 (kJ/kg)', '㶲效率 (%)', 'Location', 'northwest');
title('图3.7 系统各设备㶲损与㶲效率');
savefig('图3.7_系统各设备㶲损与㶲效率.fig');
saveas(gcf, '图3.7_系统各设备㶲损与㶲效率.png');

% 总结输出
disp('参数分析和图表生成完成！');
disp('基准工况的循环热效率为：');
[eta_base, w_net_base, ~, ~] = calculate_cycle_performance(...
    default_P1, default_P2, default_P4, default_P8, default_alpha, ...
    default_n_hp_turb, default_n_lp_turb, default_N_mc_a, default_N_mc_b, default_N_rc, ...
    default_e1, default_e2, default_T1, default_T7);
fprintf('循环热效率: %.4f%%\n', eta_base*100);
fprintf('净输出功率: %.2f kJ/kg\n', w_net_base/1000);

%% 循环性能计算函数
function [eta, W_net, exergy_destruction, exergy_efficiency] = calculate_cycle_performance(...
    P1, P2, P4, P8, alpha, n_hp_turb, n_lp_turb, N_mc_a, N_mc_b, N_rc, e1, e2, T1, T7)
    % 此函数执行超临界CO2布雷顿循环的核心计算，返回热效率、净功率和复杂度指标
    
    % 固定参数设置
    T3 = T1;               % 再热温度与初始温度相同
    P3 = P2;               % 再热器压力保持不变
    P7 = P4;               % 冷却器出口压力
    P9 = P8;               % 中间冷却压力保持不变
    T9 = T7;               % 中间冷却后温度
    P10 = P1;              % 主压缩机b出口压力
    P11 = P1;              % 副压缩机出口压力
    P6 = P4;               % 分流点压力
    P5 = P4;               % 高温回热器热侧出口压力
    
    % 环境参考状态（用于㶲分析）
    T0 = 298.15;  % 环境温度(K)
    P0 = 0.101325; % 环境压力(MPa)
    H0 = refpropm('h','T',T0,'P',P0*1000,'CO2');
    S0 = refpropm('s','T',T0,'P',P0*1000,'CO2');
    
    %------初始计算------
    % 高压透平入口
    S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
    H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');
    
    % 高压透平出口
    S2_is = S1;
    H2_is = refpropm('h','P',P2*1000,'S',S2_is,'CO2');
    H2 = H1 - (H1 - H2_is)*n_hp_turb;
    T2 = refpropm('t','P',P2*1000,'H',H2,'CO2');
    S2 = refpropm('s','T',T2,'P',P2*1000,'CO2');
    
    % 再热器出口
    H3 = refpropm('h','T',T3,'P',P3*1000,'CO2');
    S3 = refpropm('s','T',T3,'P',P3*1000,'CO2');
    
    % 低压透平出口
    S4_is = S3;
    H4_is = refpropm('h','P',P4*1000,'S',S4_is,'CO2');
    H4 = H3 - (H3 - H4_is)*n_lp_turb;
    T4 = refpropm('t','P',P4*1000,'H',H4,'CO2');
    S4 = refpropm('s','T',T4,'P',P4*1000,'CO2');
    
    % 冷却器出口/主压缩机a入口
    H7 = refpropm('h','T',T7,'P',P7*1000,'CO2');
    S7 = refpropm('s','T',T7,'P',P7*1000,'CO2');
    
    % 主压缩机a出口
    S8_is = S7;
    H8_is = refpropm('h','P',P8*1000,'S',S8_is,'CO2');
    H8 = H7 + (H8_is - H7)/N_mc_a;
    T8 = refpropm('t','P',P8*1000,'H',H8,'CO2');
    S8 = refpropm('s','T',T8,'P',P8*1000,'CO2');
    
    % 中间冷却器出口
    H9 = refpropm('h','T',T9,'P',P9*1000,'CO2');
    S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');
    
    % 主压缩机b出口
    S10_is = S9;
    H10_is = refpropm('h','P',P10*1000,'S',S10_is,'CO2');
    H10 = H9 + (H10_is - H9)/N_mc_b;
    T10 = refpropm('t','P',P10*1000,'H',H10,'CO2');
    S10 = refpropm('s','T',T10,'P',P10*1000,'CO2');
    
    %------迭代初始化------
    % 使用新方法初始化
    T5_guess = 600;
    
    % 预估副压缩机状态
    T11_est = 320;
    H11_est = refpropm('h','T',T11_est,'P',P6*1000,'CO2');
    S11_est = refpropm('s','T',T11_est,'P',P6*1000,'CO2');
    
    S12_is = S11_est;
    H12_is = refpropm('h','P',P11*1000,'S',S12_is,'CO2');
    H12 = H11_est + (H12_is - H11_est)/N_rc;
    T12 = refpropm('t','P',P11*1000,'H',H12,'CO2');
    S12 = refpropm('s','T',T12,'P',P11*1000,'CO2');
    
    % 假设合流点状态
    H13_ideal = (1-alpha)*H10 + alpha*H12;
    entropy_gen_factor = 0.03;
    entropy_loss = entropy_gen_factor * abs(H12 - H10);
    H13 = H13_ideal - entropy_loss;
    P13 = P10;
    T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
    S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
    
    % 使用效率法估计回热器温度
    P14 = P13;
    T14 = T13 + e1 * (T5_guess - T13);
    H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
    S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
    
    % 高温回热器冷侧出口温度
    P15 = P14;
    T15 = T14 + e2 * (T4 - T14);
    H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
    S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
    
    % 反向计算高温回热器热侧出口
    Q_HT_cold = H15 - H14;
    H5 = H4 - Q_HT_cold;
    T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
    S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
    
    % 反向计算低温回热器热侧出口
    Q_LT_cold = H14 - H13;
    H6 = H5 - Q_LT_cold;
    T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
    S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
    
    % 更新分流点状态
    T11 = T6;
    H11 = H6;
    S11 = S6;
    
    % 更新副压缩机出口
    S12_is = S11;
    H12_is = refpropm('h','P',P11*1000,'S',S12_is,'CO2');
    H12 = H11 + (H12_is - H11)/N_rc;
    T12 = refpropm('t','P',P11*1000,'H',H12,'CO2');
    S12 = refpropm('s','T',T12,'P',P11*1000,'CO2');
    
    % 重新计算合流点
    H13_ideal = (1-alpha)*H10 + alpha*H12;
    entropy_gen_factor = 0.03;
    entropy_loss = entropy_gen_factor * abs(H12 - H10);
    H13 = H13_ideal - entropy_loss;
    T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
    S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
    
    %------热平衡迭代求解------
    % 热平衡迭代参数
    max_iter = 50;
    relax_coef = 0.5;
    converge_tol = 0.001;
    
    % 迭代求解
    object = 1;
    iter = 0;
    
    while object && iter < max_iter
        iter = iter + 1;
        
        % 保存上一轮值
        H5_old = H5;
        H6_old = H6;
        H12_old = H12;
        H13_old = H13;
        
        % 更新高温回热器
        Q_HT_hot = H4 - H5;
        Q_HT_cold = H13 - H12;
        
        % 调整高温回热器热侧出口焓
        H5_new = H4 - Q_HT_cold;
        
        % 带松弛系数的更新
        relax = relax_coef;
        H5 = H5_old + relax*(H5_new - H5_old);
        T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
        S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
        
        % 更新低温回热器
        Q_LT_hot = H5 - H6;
        
        % 计算低温回热器冷侧出口温度
        T10_LTR_out = T10 + e1 * (T5 - T10);
        H10_LTR_out = refpropm('h','T',T10_LTR_out,'P',P10*1000,'CO2');
        S10_LTR_out = refpropm('s','T',T10_LTR_out,'P',P10*1000,'CO2');
        
        % 更新副压缩机入口状态
        T6_rc_in = T6;
        H6_rc_in = H6;
        S6_rc_in = S6;
        
        % 更新副压缩机出口状态
        S11_is = S6_rc_in;
        H11_is = refpropm('h','P',P11*1000,'S',S11_is,'CO2');
        H11 = H6_rc_in + (H11_is - H6_rc_in)/N_rc;
        T11 = refpropm('t','P',P11*1000,'H',H11,'CO2');
        S11 = refpropm('s','T',T11,'P',P11*1000,'CO2');
        
        % 计算低温回热器副路出口温度
        T11_LTR_out = T11 + e1 * (T5 - T11);
        H11_LTR_out = refpropm('h','T',T11_LTR_out,'P',P11*1000,'CO2');
        S11_LTR_out = refpropm('s','T',T11_LTR_out,'P',P11*1000,'CO2');
        
        % 总吸热量
        Q_LT_cold_main = H10_LTR_out - H10;
        Q_LT_cold_bypass = H11_LTR_out - H11;
        Q_LT_cold = (1-alpha)*Q_LT_cold_main + alpha*Q_LT_cold_bypass;
        
        % 调整低温回热器热侧出口焓
        H6_new = H5 - Q_LT_cold;
        H6 = H6_old + relax*(H6_new - H6_old);
        T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
        S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
        
        % 计算合流点
        H12_ideal = (1-alpha)*H10_LTR_out + alpha*H11_LTR_out;
        entropy_loss = entropy_gen_factor * abs(H11_LTR_out - H10_LTR_out);
        H12 = H12_ideal - entropy_loss;
        T12 = refpropm('t','P',P10*1000,'H',H12,'CO2');
        S12 = refpropm('s','T',T12,'P',P10*1000,'CO2');
        
        % 更新高温回热器冷侧出口
        T13 = T12 + e2 * (T4 - T12);
        H13 = refpropm('h','T',T13,'P',P13*1000,'CO2');
        S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
        
        % 检查收敛性
        dH5 = abs(H5 - H5_old)/abs(H5);
        dH6 = abs(H6 - H6_old)/abs(H6);
        dH12 = abs(H12 - H12_old)/abs(H12);
        dH13 = abs(H13 - H13_old)/abs(H13);
        
        if max([dH5, dH6, dH12, dH13]) < converge_tol
            object = 0;
        end
    end
    
    %------计算循环效率------
    % 功率计算
    W_hp_turb = H1 - H2;
    W_lp_turb = H3 - H4;
    W_mc_a = H8 - H7;
    W_mc_b = H10 - H9;
    W_rc = H11 - H6;
    W_comp_total = (1-alpha)*(W_mc_a + W_mc_b) + alpha*W_rc;
    W_turb_total = W_hp_turb + W_lp_turb;
    W_net = W_turb_total - W_comp_total;
    
    % 热量计算
    Q_heater = H1 - H13;
    Q_reheater = H3 - H2;
    Q_cooler = (1-alpha)*(H6 - H7);
    Q_intercooler = (1-alpha)*(H8 - H9);
    Q_in_total = Q_heater + Q_reheater;
    
    % 循环热效率
    eta = W_net / Q_in_total;
    
    %------㶲分析计算------
    % 计算各点㶲值 (kJ/kg)
    E1 = (H1 - H0)/1000 - T0*(S1 - S0)/1000;
    E2 = (H2 - H0)/1000 - T0*(S2 - S0)/1000;
    E3 = (H3 - H0)/1000 - T0*(S3 - S0)/1000;
    E4 = (H4 - H0)/1000 - T0*(S4 - S0)/1000;
    E5 = (H5 - H0)/1000 - T0*(S5 - S0)/1000;
    E6 = (H6 - H0)/1000 - T0*(S6 - S0)/1000;
    E7 = (H7 - H0)/1000 - T0*(S7 - S0)/1000;
    E8 = (H8 - H0)/1000 - T0*(S8 - S0)/1000;
    E9 = (H9 - H0)/1000 - T0*(S9 - S0)/1000;
    E10 = (H10 - H0)/1000 - T0*(S10 - S0)/1000;
    E11 = (H11 - H0)/1000 - T0*(S11 - S0)/1000;
    E12 = (H12 - H0)/1000 - T0*(S12 - S0)/1000;
    E13 = (H13 - H0)/1000 - T0*(S13 - S0)/1000;
    
    % 计算等熵过程㶲值 (kJ/kg)
    E2_is = (H2_is - H0)/1000 - T0*(S1 - S0)/1000;
    E4_is = (H4_is - H0)/1000 - T0*(S3 - S0)/1000;
    E8_is = (H8_is - H0)/1000 - T0*(S7 - S0)/1000;
    E10_is = (H10_is - H0)/1000 - T0*(S9 - S0)/1000;
    E11_is = (H11_is - H0)/1000 - T0*(S11 - S0)/1000;
    
    % 计算㶲损
    % 高压透平㶲损
    Ed_hp_turb = E1 - E2 - W_hp_turb/1000;
    % 低压透平㶲损
    Ed_lp_turb = E3 - E4 - W_lp_turb/1000;
    % 主压缩机a㶲损
    Ed_mc_a = (1-alpha)*(W_mc_a/1000 - (E8 - E7));
    % 主压缩机b㶲损
    Ed_mc_b = (1-alpha)*(W_mc_b/1000 - (E10 - E9));
    % 副压缩机㶲损
    Ed_rc = alpha*(W_rc/1000 - (E11 - E6));
    % 高温回热器㶲损
    Ed_ht_recup = (1-alpha)*(E4 - E5) - (E13 - E12);
    % 低温回热器㶲损
    Ed_lt_recup = (1-alpha)*(E5 - E6) - (E12 - E10);
    % 冷却器㶲损
    Ed_cooler = (1-alpha)*E6;
    % 合流点㶲损
    Ed_mixing = (1-alpha)*E10 + alpha*E11 - E12;
    % 加热器㶲损（简化计算）
    Ed_heater = Q_heater/1000 * (1 - T0/T1) - (E1 - E13);
    
    % 各设备㶲损合集
    exergy_destruction = [Ed_hp_turb, Ed_lp_turb, Ed_mc_a, Ed_mc_b, Ed_rc, ...
                         Ed_ht_recup, Ed_lt_recup, Ed_cooler, Ed_mixing, Ed_heater];
    
    % 计算各设备㶲效率
    % 高压透平㶲效率
    Ee_hp_turb = (W_hp_turb/1000) / (E1 - E2_is);
    % 低压透平㶲效率
    Ee_lp_turb = (W_lp_turb/1000) / (E3 - E4_is);
    % 主压缩机a㶲效率
    Ee_mc_a = (E8_is - E7) / (W_mc_a/1000);
    % 主压缩机b㶲效率
    Ee_mc_b = (E10_is - E9) / (W_mc_b/1000);
    % 副压缩机㶲效率
    Ee_rc = (E11_is - E6) / (W_rc/1000);
    % 高温回热器㶲效率
    Ee_ht_recup = (E13 - E12) / ((1-alpha)*(E4 - E5));
    % 低温回热器㶲效率
    Ee_lt_recup = (E12 - E10) / ((1-alpha)*(E5 - E6));
    % 冷却器㶲效率（通常很低）
    Ee_cooler = 0;
    % 合流点㶲效率
    Ee_mixing = E12 / ((1-alpha)*E10 + alpha*E11);
    % 加热器㶲效率
    Ee_heater = (E1 - E13) / (Q_heater/1000 * (1 - T0/T1));
    
    % 各设备㶲效率合集
    exergy_efficiency = [Ee_hp_turb, Ee_lp_turb, Ee_mc_a, Ee_mc_b, Ee_rc, ...
                         Ee_ht_recup, Ee_lt_recup, Ee_cooler, Ee_mixing, Ee_heater] * 100; % 转为百分比
end 