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

% 安全常数定义
T_CO2_min = 220; % CO₂最低安全温度(K)
min_temp_diff = 10; % 最小温差(K)
max_pressure = 30; % 最大安全压力(MPa)

% 安全检查：压力参数
P_high = min(P_high, max_pressure);
P_low = min(P_low, P_high - 1);
P_reheat = min(P_reheat, P_high);
P_reheat = max(P_reheat, P_low + 1);
P_intercool = min(P_intercool, P_high - 1);
P_intercool = max(P_intercool, P_low + 1);

% 创建状态点结构体数组并初始化所有字段为0
for i = 1:17
    states(i) = struct('T', 0, 'P', 0, 'h', 0, 's', 0, 'rho', 0, 'cp', 0);
end

% 初始化关键状态点
% 状态点1 (高压透平入口)
states(1).T = T_high;
states(1).P = P_high;
states(1).h = refpropm('H', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000;
states(1).s = refpropm('S', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000;
states(1).rho = refpropm('D', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2');
states(1).cp = refpropm('C', 'T', states(1).T, 'P', states(1).P * 1000, 'CO2') / 1000;

% 状态点8 (冷却器出口/主压缩机a入口)
states(8).T = T_low;
states(8).P = P_low;
states(8).h = refpropm('H', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;
states(8).s = refpropm('S', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;
states(8).rho = refpropm('D', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2');
states(8).cp = refpropm('C', 'T', states(8).T, 'P', states(8).P * 1000, 'CO2') / 1000;

% 初始估计状态点14 (合流点出口/低温回热器冷侧入口)
states(14).T = max(T_low + 50, T_CO2_min); % 确保温度高于CO₂的最低允许温度
states(14).P = P_high;
states(14).h = refpropm('H', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
states(14).s = refpropm('S', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
states(14).rho = refpropm('D', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2');
states(14).cp = refpropm('C', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;

% 计算各状态点热力学参数的最大迭代次数和收敛容差
max_iter = 100;
tol = 1e-6;

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
    % 安全检查：确保温度不低于最低允许温度
    states(2).T = max(states(2).T, T_CO2_min);
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
    % 安全检查：确保温度不低于最低允许温度
    states(4).T = max(states(4).T, T_CO2_min);
    states(4).s = refpropm('S', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    states(4).rho = refpropm('D', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2');      % kg/m³
    states(4).cp = refpropm('C', 'T', states(4).T, 'P', states(4).P * 1000, 'CO2') / 1000; % kJ/(kg·K)
    
    % ========== 2. 回热器计算（修改后的逻辑） ==========
    % ---------- 2.1 高温回热器(HT)计算 ----------
    % a. 确定热侧进出口状态
    % 状态点4已经计算好（高温回热器热侧入口）
    
    % ---- HT 回热器调试输出 ----
    Q_HT_hot      = m_dot * (states(4).h - states(5).h) * eta_recup_HT;
    T16_ideal     = states(4).T - HT_dT;
    T16_minByDiff = states(15).T + min_temp_diff;
    h16_ideal     = refpropm('H','T', max(T16_ideal, T16_minByDiff), 'P', states(15).P*1000, 'CO2')/1000;
    Q_HT_cold_max = m_dot * (h16_ideal - states(15).h);
    Q_HT_actual   = min(Q_HT_hot, Q_HT_cold_max);
    fprintf('HT recup iter %d: h4=%.1f, h5=%.1f, Q_hot=%.1f, Q_cold_max=%.1f, Q_act=%.1f, T16_i=%.1f, T15=%.1f\n', ...
            iter, states(4).h, states(5).h, Q_HT_hot, Q_HT_cold_max, Q_HT_actual, T16_ideal, states(15).T);
    % -----------------------------
    
    % 使用端差约束确定热侧出口温度，但要确保温度不低于最低安全温度
    % 假设冷侧入口温度states(15).T（暂不知道，使用上次迭代结果或初始估计）
    if iter == 1
        % 第一次迭代，使用初始估计
        states(15).T = max(states(4).T - 100, T_CO2_min);
        states(15).P = P_high;
        states(15).h = refpropm('H', 'T', states(15).T, 'P', states(15).P * 1000, 'CO2') / 1000;
    end
    
    % 确定热侧出口温度(状态点5)，确保有足够的温差
    states(5).T = max(states(15).T + min_temp_diff, T_CO2_min);
    states(5).P = states(4).P; % 忽略压降
    states(5).h = refpropm('H', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2') / 1000;
    states(5).s = refpropm('S', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2') / 1000;
    states(5).rho = refpropm('D', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2');
    states(5).cp = refpropm('C', 'T', states(5).T, 'P', states(5).P * 1000, 'CO2') / 1000;
    
    % b. 计算高温回热器热侧热负荷
    Q_HT_available = Q_HT_hot * eta_recup_HT;
    
    % c. 设定冷侧出口温度，需满足端差约束
    % 理想情况下，冷侧出口温度 = 热侧入口温度 - 端差
    T16_ideal = states(4).T - HT_dT;
    % 考虑安全限制
    T16_ideal = max(T16_ideal, states(15).T + min_temp_diff);
    T16_ideal = max(T16_ideal, T_CO2_min);
    
    % 计算理想情况下的冷侧出口状态
    states(16).P = states(15).P; % 忽略压降
    states(16).T = T16_ideal;
    % 安全检查：确保焓值在合理范围内
    try
        h16_ideal = refpropm('H', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    catch
        % 如果计算出错，使用安全值
        states(16).T = max(states(15).T + min_temp_diff, T_CO2_min);
        h16_ideal = refpropm('H', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    end
    
    % d. 计算冷侧实际吸收热量
    Q_HT_cold_max = m_dot * (h16_ideal - states(15).h);
    % 取热侧释放热量与冷侧最大吸收热量的较小值
    Q_HT_actual = min(Q_HT_available, Q_HT_cold_max);
    
    % e. 根据实际换热量计算状态点16的最终状态
    states(16).h = states(15).h + Q_HT_actual / m_dot;
    % 安全检查：确保焓值不低于最小值
    try
        min_h16 = refpropm('H', 'T', T_CO2_min, 'P', states(16).P * 1000, 'CO2') / 1000;
        states(16).h = max(states(16).h, min_h16);
        states(16).T = refpropm('T', 'P', states(16).P * 1000, 'H', states(16).h * 1000, 'CO2');
    catch
        % 如果计算出错，保持温度不变，重新计算焓值
        states(16).T = max(states(15).T + min_temp_diff, T_CO2_min);
        states(16).h = refpropm('H', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    end
    states(16).s = refpropm('S', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    states(16).rho = refpropm('D', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2');
    states(16).cp = refpropm('C', 'T', states(16).T, 'P', states(16).P * 1000, 'CO2') / 1000;
    
    % ---------- 2.2 低温回热器(LT)计算（改进） ----------
    % 目标：同时满足能量平衡 Q_hot = Q_cold 和端差约束 T5 - T6 = LT_dT
    % 定义待求解变量 T5（热侧出口温度），T6 = T5 - LT_dT
    fun = @(T5) m_dot*(refpropm('H','T',T5,'P',states(5).P*1000,'CO2')/1000 ...
                     - refpropm('H','T',T5-LT_dT,'P',states(5).P*1000,'CO2')/1000)*eta_recup_LT ...
               - m_dot*(refpropm('H','T',T5-LT_dT,'P',states(5).P*1000,'CO2')/1000 - states(14).h);
    T5_guess = states(14).T + LT_dT*1.5;
    % 使用 fzero 求解 T5
    T5 = fzero(fun, T5_guess);
    T6 = T5 - LT_dT;
    % 更新状态点5和6
    states(5).T = T5; states(6).T = T6;
    % 同步压力，忽略压降
    states(6).P = states(5).P;
    states(5).h = refpropm('H','T',T5,'P',states(5).P*1000,'CO2')/1000;
    states(6).h = refpropm('H','T',T6,'P',states(6).P*1000,'CO2')/1000;
    states(5).s = refpropm('S','T',T5,'P',states(5).P*1000,'CO2')/1000;
    states(6).s = refpropm('S','T',T6,'P',states(6).P*1000,'CO2')/1000;
    states(5).rho = refpropm('D','T',T5,'P',states(5).P*1000,'CO2');
    states(6).rho = refpropm('D','T',T6,'P',states(6).P*1000,'CO2');
    states(5).cp = refpropm('C','T',T5,'P',states(5).P*1000,'CO2')/1000;
    states(6).cp = refpropm('C','T',T6,'P',states(6).P*1000,'CO2')/1000;
    % 计算实际低温回热器热负荷
    Q_recup_LT = m_dot*(states(5).h - states(6).h)*eta_recup_LT;

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
    % 安全检查：确保温度不低于最低允许温度
    states(9).T = max(states(9).T, T_CO2_min);
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
    % 安全检查：确保温度不低于最低允许温度
    states(11).T = max(states(11).T, T_CO2_min);
    states(11).s = refpropm('S', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2') / 1000;
    states(11).rho = refpropm('D', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2');
    states(11).cp = refpropm('C', 'T', states(11).T, 'P', states(11).P * 1000, 'CO2') / 1000;
    
    % 副路 - 副压缩机
    % 状态点13 (副压缩机出口/合流点)
    % ---------- 副压缩机理想等熵出口（限幅处理） ----------
    % 查询该压力下熵的上下限（kJ/kg·K）
    Smin = refpropm('S','T',T_low,'P',P_high*1000,'CO2')/1000;
    Smax = refpropm('S','T',T_high,'P',P_high*1000,'CO2')/1000;
    % 原始熵值并限幅到 [Smin,Smax]
    s13s_unclamped = states(12).s * 1000;
    s13s_clamped   = min(max(s13s_unclamped, Smin*1000), Smax*1000);
    % 安全调用 REFPROP 计算焓，失败时降级使用等温焓
    try
        h13s = refpropm('H','P',P_high*1000,'S', s13s_clamped,'CO2')/1000;
    catch
        warning('副压缩机等熵变换熵超范围，降级使用等温焓');
        h13s = refpropm('H','T', states(12).T, 'P', P_high*1000,'CO2')/1000;
    end
    % 实际焓变 = 等熵焓变 / 压缩机效率
    states(13).h = states(12).h + (h13s - states(12).h) / eta_c_recomp;
    states(13).P = P_high; % MPa
    states(13).T = refpropm('T', 'P', states(13).P * 1000, 'H', states(13).h * 1000, 'CO2');
    % 安全检查：确保温度不低于最低允许温度
    states(13).T = max(states(13).T, T_CO2_min);
    states(13).s = refpropm('S', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2') / 1000;
    states(13).rho = refpropm('D', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2');
    states(13).cp = refpropm('C', 'T', states(13).T, 'P', states(13).P * 1000, 'CO2') / 1000;
    
    % 合流 - 状态点14 (合流点出口/低温回热器冷侧入口)
    states(14).h = (1 - split_ratio) * states(11).h + split_ratio * states(13).h;
    states(14).P = P_high; % MPa
    % 使用焓值和压力计算温度，并确保在有效范围内
    try
        states(14).T = refpropm('T', 'P', states(14).P * 1000, 'H', states(14).h * 1000, 'CO2');
        % 安全检查：确保温度不低于最低允许温度
        states(14).T = max(states(14).T, T_CO2_min);
    catch
        % 如果计算出错，使用安全值
        states(14).T = max((states(11).T + states(13).T) / 2, T_CO2_min);
        states(14).h = refpropm('H', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    end
    states(14).s = refpropm('S', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    states(14).rho = refpropm('D', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2');
    states(14).cp = refpropm('C', 'T', states(14).T, 'P', states(14).P * 1000, 'CO2') / 1000;
    
    % ========== 5. 加热器计算 ==========
    % 状态点17 (加热器出口/高压透平入口) - 与状态点1相同
    states(17) = states(1);
    
    % 检查收敛性 - 比较重要状态点的温度
    diff_max = 0;
    key_points = [1, 4, 6, 14, 15, 16]; % 选择关键状态点进行收敛性检查
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
Q_recup_HT = m_dot * (states(16).h - states(15).h); % 高温回热器热负荷，kW
Q_recup_LT = m_dot * (states(15).h - states(14).h); % 低温回热器热负荷，kW

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
performance.Q_recup_HT = Q_recup_HT;
performance.Q_recup_LT = Q_recup_LT;
performance.W_turbine = W_turbine;
performance.W_compressor = W_compressor;
performance.W_net = W_net;
performance.Q_in = Q_in;
performance.eta_thermal = eta_thermal;
performance.mass_flow = m_dot;
performance.split_ratio = split_ratio;
performance.iter_count = iter;

% 迭代结束，打印全局能量流
fprintf('\n===== 全局能量流 (kW) =====\n');
fprintf('Q_heater   = %.1f\n', m_dot*(states(1).h - states(16).h));
fprintf('Q_reheater = %.1f\n', m_dot*(states(3).h - states(2).h));
fprintf('W_turbine  = %.1f\n', performance.W_turbine);
fprintf('W_comp     = %.1f\n', performance.W_compressor);
fprintf('Q_cooler   = %.1f\n', performance.Q_cooler);
fprintf('Q_interc   = %.1f\n', performance.Q_intercooler);
fprintf('Balance    = Q_in - (W_net + Q_cooler + Q_interc) = %.1f kW\n\n', ...
        performance.Q_in - (performance.W_net + performance.Q_cooler + performance.Q_intercooler));

end