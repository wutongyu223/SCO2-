function [states, performance] = calculate_cycle(params)
% CALCULATE_CYCLE 计算超临界CO₂布雷顿循环的所有状态点和性能指标
%
% 输入参数:
%   params - 包含以下字段的结构体:
%     - P_high: 最高压力 (MPa)
%     - P_low: 最低压力 (MPa)
%     - T_high: 最高温度 (K)
%     - T_low: 最低温度 (K)
%     - P_reheat: 再热压力 (MPa)
%     - P_intercool: 中间冷却压力 (MPa)
%     - split_ratio: 分流比例
%     - HT_recuperator_dT: 高温回热器端差 (K)
%     - LT_recuperator_dT: 低温回热器端差 (K)
%     - eta_turbine_HP: 高压透平效率
%     - eta_turbine_LP: 低压透平效率
%     - eta_compressor_main: 主压缩机效率
%     - eta_compressor_recomp: 副压缩机效率
%     - eta_recuperator_HT: 高温回热器效率
%     - eta_recuperator_LT: 低温回热器效率
%     - eta_heater: 加热器效率
%     - m_dot: 工质质量流量 (kg/s)
%
% 输出参数:
%   states - 包含所有状态点热力学参数的结构体数组(1-17)
%     每个状态点包含: T, P, h, s, rho, cp 等参数
%   performance - 循环性能指标的结构体
%     包含: 热效率, 净功率, 各组件功率/热负荷等

% 创建状态点结构体数组
states = struct('T',zeros(1,17),'P',zeros(1,17),'h',zeros(1,17),...
                's',zeros(1,17),'rho',zeros(1,17),'cp',zeros(1,17));

% 提取参数
P_high = params.P_high;           % MPa
P_low = params.P_low;             % MPa
T_high = params.T_high;           % K
T_low = params.T_low;             % K
P_reheat = params.P_reheat;       % MPa
P_intercool = params.P_intercool; % MPa
split_ratio = params.split_ratio; % 分流比例
HT_dT = params.HT_recuperator_dT; % K
LT_dT = params.LT_recuperator_dT; % K
eta_t_HP = params.eta_turbine_HP;
eta_t_LP = params.eta_turbine_LP;
eta_c_main = params.eta_compressor_main;
eta_c_recomp = params.eta_compressor_recomp;
eta_recup_HT = params.eta_recuperator_HT;
eta_recup_LT = params.eta_recuperator_LT;
eta_heater = params.eta_heater;
m_dot = params.m_dot;             % kg/s

% 初始化状态点1 (高压透平入口)
states(1).T = T_high;
states(1).P = P_high;

% 初始化状态点7 (分流后主路冷却器入口) - 假设先用T_low + LT_dT估计
states(7).T = T_low + LT_dT; 
states(7).P = P_low;

% 初始化状态点8 (冷却器出口/主压缩机a入口)
states(8).T = T_low;
states(8).P = P_low;

% 初始化状态点12 (分流后副路副压缩机入口) - 与状态点6相同，但先初始化为与状态点7相同
states(12).T = states(7).T;
states(12).P = states(7).P;

% 计算各状态点热力学参数的最大迭代次数和收敛容差
max_iter = 100;
tol = 1e-4;

% 主循环计算
for iter = 1:max_iter
    % 保存上一次迭代的状态点，用于检查收敛性
    states_prev = states;
    
    % ========== 1. 透平和再热器计算 ==========
    % 状态点1 (高压透平入口)
    states(1).h = refpropm('H', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000; % kJ/kg
    states(1).s = refpropm('S', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    states(1).rho = refpropm('D', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2');      % kg/m³
    states(1).cp = refpropm('C', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    
    % 状态点2 (高压透平出口/再热器入口) - 等熵膨胀后使用透平效率计算
    s2s = states(1).s; % 等熵过程
    h2s = refpropm('H', 'P', P_reheat * 1000, 'S', s2s * 1000, 'CO2') / 1000; % kJ/kg
    % 实际焓变 = 等熵焓变 / 透平效率
    states(2).h = states(1).h - eta_t_HP * (states(1).h - h2s); % kJ/kg
    states(2).P = P_reheat; % MPa
    states(2).T = refpropm('T', 'P', states(2).P * 1000, 'H', states(2).h * 1000, 'CO2'); % K
    states(2).s = refpropm('S', 'T', states(2).T, 'P', states(2).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    states(2).rho = refpropm('D', 'T', states(2).T, 'P', states(2).P * 1000, 'CO2');      % kg/m³
    states(2).cp = refpropm('C', 'T', states(2).T, 'P', states(2).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    
    % 状态点3 (再热器出口/低压透平入口)
    states(3).T = T_high; % 再热到最高温度
    states(3).P = states(2).P; % 压力不变
    states(3).h = refpropm('H', 'T', states(3).T, 'P', states(3).P * 1000, 'CO2') / 1000; % kJ/kg
    states(3).s = refpropm('S', 'T', states(3).T, 'P', states(3).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    states(3).rho = refpropm('D', 'T', states(3).T, 'P', states(3).P * 1000, 'CO2');      % kg/m³
    states(3).cp = refpropm('C', 'T', states(3).T, 'P', states(3).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    
    % 状态点4 (低压透平出口/高温回热器热侧入口)
    s4s = states(3).s; % 等熵过程
    h4s = refpropm('H', 'P', P_low * 1000, 'S', s4s * 1000, 'CO2') / 1000; % kJ/kg
    % 实际焓变 = 等熵焓变 / 透平效率
    states(4).h = states(3).h - eta_t_LP * (states(3).h - h4s); % kJ/kg
    states(4).P = P_low; % MPa
    states(4).T = refpropm('T', 'P', states(4).P * 1000, 'H', states(4).h * 1000, 'CO2'); % K
    states(4).s = refpropm('S', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    states(4).rho = refpropm('D', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2');      % kg/m³
    states(4).cp = refpropm('C', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    
    % ========== 2. 回热器计算 (需要迭代) ==========
    % 这里我们先假设知道状态点16 (从上一轮迭代)，计算回热器热侧
    
    % 高温回热器热侧计算
    % 状态点5 (高温回热器热侧出口/低温回热器热侧入口)
    if iter == 1
        % 第一次迭代时，使用端差估计状态点5的温度
        states(16).T = states(4).T - HT_dT;
        states(16).P = P_high;
        states(16).h = refpropm('H', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    end
    
    % 高温回热器热负荷
    Q_recup_HT = eta_recup_HT * m_dot * (states(16).h - states(15).h);
    states(5).h = states(4).h - Q_recup_HT / m_dot;
    states(5).P = states(4).P; % 忽略压降
    states(5).T = refpropm('T', 'P', states(5).P * 1000, 'H', states(5).h * 1000, 'CO2');
    states(5).s = refpropm('S', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2') / 1000;
    states(5).rho = refpropm('D', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2');
    states(5).cp = refpropm('C', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2') / 1000;
    
    % 低温回热器热侧计算
    % 状态点6 (低温回热器热侧出口/分流点)
    if iter == 1
        % 第一次迭代时，使用端差估计状态点14的温度
        states(14).T = states(7).T + LT_dT;
        states(14).P = P_intercool; % 假设是中间冷却压力
        states(14).h = refpropm('H', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    end
    
    % 低温回热器热负荷
    Q_recup_LT = eta_recup_LT * m_dot * (states(15).h - states(14).h);
    states(6).h = states(5).h - Q_recup_LT / m_dot;
    states(6).P = states(5).P; % 忽略压降
    states(6).T = refpropm('T', 'P', states(6).P * 1000, 'H', states(6).h * 1000, 'CO2');
    states(6).s = refpropm('S', 'T', states(6).T, 'P', states(6).P * 1000, 'CO2') / 1000;
    states(6).rho = refpropm('D', 'T', states(6).T, 'P', states(6).P * 1000, 'CO2');
    states(6).cp = refpropm('C', 'T', states(6).T, 'P', states(6).P * 1000, 'CO2') / 1000;
    
    % ========== 3. 分流、冷却器和压缩机计算 ==========
    % 分流 - 状态点7和状态点12的热力学性质与状态点6相同，仅流量不同
    states(7) = states(6);
    states(12) = states(6);
    
    % 主路 - 冷却器和主压缩机
    % 状态点8 (冷却器出口/主压缩机a入口)
    states(8).T = T_low;
    states(8).P = P_low;
    states(8).h = refpropm('H', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;
    states(8).s = refpropm('S', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;
    states(8).rho = refpropm('D', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2');
    states(8).cp = refpropm('C', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;
    
    % 状态点9 (主压缩机a出口/中间冷却入口)
    s9s = states(8).s; % 等熵过程
    h9s = refpropm('H', 'P', P_intercool * 1000, 'S', s9s * 1000, 'CO2') / 1000;
    % 实际焓变 = 等熵焓变 / 压缩机效率
    states(9).h = states(8).h + (h9s - states(8).h) / eta_c_main;
    states(9).P = P_intercool; % MPa
    states(9).T = refpropm('T', 'P', states(9).P * 1000, 'H', states(9).h * 1000, 'CO2');
    states(9).s = refpropm('S', 'T', states(9).T, 'P', states(9).P * 1000, 'CO2') / 1000;
    states(9).rho = refpropm('D', 'T', states(9).T, 'P', states(9).P * 1000, 'CO2');
    states(9).cp = refpropm('C', 'T', states(9).T, 'P', states(9).P * 1000, 'CO2') / 1000;
    
    % 状态点10 (中间冷却出口/主压缩机b入口)
    states(10).T = T_low;
    states(10).P = states(9).P; % 忽略压降
    states(10).h = refpropm('H', 'T', states(10).T, 'P', states(10).P * 1000, 'CO2') / 1000;
    states(10).s = refpropm('S', 'T', states(10).T, 'P', states(10).P * 1000, 'CO2') / 1000;
    states(10).rho = refpropm('D', 'T', states(10).T, 'P', states(10).P * 1000, 'CO2');
    states(10).cp = refpropm('C', 'T', states(10).T, 'P', states(10).P * 1000, 'CO2') / 1000;
    
    % 状态点11 (主压缩机b出口/合流点)
    s11s = states(10).s; % 等熵过程
    h11s = refpropm('H', 'P', P_high * 1000, 'S', s11s * 1000, 'CO2') / 1000;
    % 实际焓变 = 等熵焓变 / 压缩机效率
    states(11).h = states(10).h + (h11s - states(10).h) / eta_c_main;
    states(11).P = P_high; % MPa
    states(11).T = refpropm('T', 'P', states(11).P * 1000, 'H', states(11).h * 1000, 'CO2');
    states(11).s = refpropm('S', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2') / 1000;
    states(11).rho = refpropm('D', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2');
    states(11).cp = refpropm('C', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2') / 1000;
    
    % 副路 - 副压缩机
    % 状态点13 (副压缩机出口/合流点)
    s13s = states(12).s; % 等熵过程
    h13s = refpropm('H', 'P', P_high * 1000, 'S', s13s * 1000, 'CO2') / 1000;
    % 实际焓变 = 等熵焓变 / 压缩机效率
    states(13).h = states(12).h + (h13s - states(12).h) / eta_c_recomp;
    states(13).P = P_high; % MPa
    states(13).T = refpropm('T', 'P', states(13).P * 1000, 'H', states(13).h * 1000, 'CO2');
    states(13).s = refpropm('S', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2') / 1000;
    states(13).rho = refpropm('D', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2');
    states(13).cp = refpropm('C', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2') / 1000;
    
    % 合流 - 状态点14 (合流点出口/低温回热器冷侧入口)
    % 能量平衡: m_14 * h_14 = m_11 * h_11 + m_13 * h_13
    states(14).h = (1 - split_ratio) * states(11).h + split_ratio * states(13).h;
    states(14).P = P_high; % MPa
    states(14).T = refpropm('T', 'P', states(14).P * 1000, 'H', states(14).h * 1000, 'CO2');
    states(14).s = refpropm('S', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    states(14).rho = refpropm('D', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2');
    states(14).cp = refpropm('C', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    
    % ========== 4. 回热器冷侧计算 ==========
    % 状态点15 (低温回热器冷侧出口/高温回热器冷侧入口)
    states(15).h = states(14).h + Q_recup_LT / m_dot;
    states(15).P = states(14).P; % 忽略压降
    states(15).T = refpropm('T', 'P', states(15).P * 1000, 'H', states(15).h * 1000, 'CO2');
    states(15).s = refpropm('S', 'T', states(15).T, 'P', states(15).P * 1000, 'CO2') / 1000;
    states(15).rho = refpropm('D', 'T', states(15).T, 'P', states(15).P * 1000, 'CO2');
    states(15).cp = refpropm('C', 'T', states(15).T, 'P', states(15).P * 1000, 'CO2') / 1000;
    
    % 状态点16 (高温回热器冷侧出口/加热器入口)
    states(16).h = states(15).h + Q_recup_HT / m_dot;
    states(16).P = states(15).P; % 忽略压降
    states(16).T = refpropm('T', 'P', states(16).P * 1000, 'H', states(16).h * 1000, 'CO2');
    states(16).s = refpropm('S', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    states(16).rho = refpropm('D', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2');
    states(16).cp = refpropm('C', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    
    % ========== 5. 加热器计算 ==========
    % 状态点17 (加热器出口/高压透平入口) - 与状态点1相同
    states(17) = states(1);
    
    % 检查收敛性 - 比较重要状态点的温度
    diff_max = 0;
    key_points = [1, 4, 6, 14, 16]; % 选择关键状态点进行收敛性检查
    for i = key_points
        if iter > 1
            diff = abs(states(i).T - states_prev(i).T);
            diff_max = max(diff_max, diff);
        end
    end
    
    if iter > 1 && diff_max < tol
        break; % 收敛，退出迭代
    end
end

% 计算系统性能
% 各组件功率/热负荷
W_turbine_HP = m_dot * (states(1).h - states(2).h); % 高压透平功率，kW
W_turbine_LP = m_dot * (states(3).h - states(4).h); % 低压透平功率，kW
W_compressor_main_a = (1 - split_ratio) * m_dot * (states(9).h - states(8).h); % 主压缩机a功耗，kW
W_compressor_main_b = (1 - split_ratio) * m_dot * (states(11).h - states(10).h); % 主压缩机b功耗，kW
W_compressor_recomp = split_ratio * m_dot * (states(13).h - states(12).h); % 副压缩机功耗，kW
Q_cooler = (1 - split_ratio) * m_dot * (states(7).h - states(8).h); % 冷却器热负荷，kW
Q_intercooler = (1 - split_ratio) * m_dot * (states(9).h - states(10).h); % 中间冷却器热负荷，kW
Q_heater = m_dot * (states(1).h - states(16).h); % 加热器热负荷，kW
Q_reheater = m_dot * (states(3).h - states(2).h); % 再热器热负荷，kW

% 系统净功率和热效率
W_turbine = W_turbine_HP + W_turbine_LP; % 总透平功率，kW
W_compressor = W_compressor_main_a + W_compressor_main_b + W_compressor_recomp; % 总压缩机功耗，kW
W_net = W_turbine - W_compressor; % 净功率，kW
Q_in = Q_heater + Q_reheater; % 总热输入，kW
eta_thermal = W_net / Q_in; % 热效率

% 输出性能指标
performance = struct();
performance.W_turbine_HP = W_turbine_HP;
performance.W_turbine_LP = W_turbine_LP;
performance.W_compressor_main_a = W_compressor_main_a;
performance.W_compressor_main_b = W_compressor_main_b;
performance.W_compressor_recomp = W_compressor_recomp;
performance.Q_cooler = Q_cooler;
performance.Q_intercooler = Q_intercooler;
performance.Q_heater = Q_heater;
performance.Q_reheater = Q_reheater;
performance.W_turbine = W_turbine;
performance.W_compressor = W_compressor;
performance.W_net = W_net;
performance.Q_in = Q_in;
performance.eta_thermal = eta_thermal;
performance.mass_flow = m_dot;
performance.split_ratio = split_ratio;
performance.iter_count = iter;

end