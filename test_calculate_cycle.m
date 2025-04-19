% test_calculate_cycle.m
clear; clc;

% 1. 构造 para
para.P_high       = 16e6;
para.P_low        = 7.5e6;
para.T_high       = 840;
para.T_low        = 305;
para.eta_t_HP     = 0.93;
para.eta_t_LP     = 0.93;
para.eta_c_main   = 0.89;
para.eta_c_recomp = 0.89;
para.eta_recup_HT = 0.86;
para.eta_recup_LT = 0.86;
para.eta_heater   = 0.94;
para.eta_reheater = 0.94; % 再热器热效率
para.m_dot        = 100;
para.P_reheat     = 10e6;
para.P_intercool  = 10e6;
para.alpha        = 0.3;
para.deltaT_HT    = 10;
para.deltaT_LT    = 10;

% 2. 构造参考态 REF
REF.T0 = 298.15;
REF.P0 = 101.325;  % kPa
REF.h0 = refpropm('H','T',REF.T0,'P',REF.P0,'CO2');
REF.s0 = refpropm('S','T',REF.T0,'P',REF.P0,'CO2');

% 3. 调用 calculate_cycle，获取状态和性能
[state, perf] = calculate_cycle(para);

% 4. 计算能量平衡误差: Q_in = W_net + Q_cool
energy_in  = perf.Q_in;
energy_out = perf.W_net + perf.Q_cool;
energy_error = abs((energy_in - energy_out) / energy_in);

% 5. 打印结果
fprintf('----------------- calculate_cycle 调试结果 -----------------\n');
fprintf('迭代状态: status = %d (0=OK,1=propertyFail,2=diverge)\n', perf.status);
fprintf('异常信息: %s\n', perf.errorMsg);
fprintf('热输入 Q_in   = %.3f kW\n', perf.Q_in);
fprintf('冷却负荷 Q_cool = %.3f kW\n', perf.Q_cool);
fprintf('净功率 W_net = %.3f kW\n', perf.W_net);
fprintf('热效率 eta_th = %.4f\n', perf.eta_th);
fprintf('能量平衡误差 Energy error = %.3f%%\n', energy_error * 100); 