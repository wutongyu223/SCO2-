%% 定义问题数据结构
% 将计算脚本的路径添加到数据结构中
data.scriptPath = 'calautatlion.m'; 
data.P_hp_in = 21; % 高压透平入口压力，MPa
data.P_lp_out = 7.5729; % 低压透平出口压力，MPa

%% 定义评价函数
function [x_repaired, objectives, constraints] = evalSCO2Cycle(x, data)
    % 修复决策变量到合理范围内
    x_repaired = x;
    
    % 解码决策变量
    P_reheat = x(1);     % 再热压力，MPa
    P_intercool = x(2);  % 中间冷却压力，MPa
    alpha = x(3);        % 分流比例
    dT_RHX_HT = x(4);    % 高温回热器端差，K
    dT_RHX_LT = x(5);    % 低温回热器端差，K
    
    % 保存当前工作区中的原始值
    try
        originalP2 = evalin('base', 'P2');
        originalP9 = evalin('base', 'P9');
        originalAlpha = evalin('base', 'alpha');
        originalTHTR = evalin('base', 'THTR');
        originalTLTR = evalin('base', 'TLTR');
    catch
        % 如果变量不存在，使用默认值
        originalP2 = 15;
        originalP9 = 12;
        originalAlpha = 0.3333;
        originalTHTR = 20;
        originalTLTR = 20;
    end
    
    % 设置新参数到工作区
    assignin('base', 'P2', P_reheat);
    assignin('base', 'P9', P_intercool);
    assignin('base', 'alpha', alpha);
    assignin('base', 'THTR', dT_RHX_HT);
    assignin('base', 'TLTR', dT_RHX_LT);
    
    try
        % 运行计算脚本并抑制命令窗口输出
        evalc('run(data.scriptPath)');
        
        % 获取计算结果
        eta = evalin('base', 'eta');
        W_net = evalin('base', 'W_net');
        
        % 简化的成本模型
        C_ref = 1e6; % USD
        Size_ref = 100e6; % W (100 MW)
        Size = W_net * 1000; % 扩大到合适的规模
        
        % 计算设备成本（使用规模修正）
        C_0 = C_ref * (Size/Size_ref)^0.6;
        
        % 计算单位成本（USD/MW）
        cost_per_MW = C_0 / (Size/1e6);
        
        % 目标函数（两个目标都是最小化）
        objectives = [-eta*100, cost_per_MW/1e6]; % 效率为百分比，成本为百万USD/MW
        
        % 约束函数
        constraints = [
            P_reheat - data.P_hp_in,        % 再热压力必须小于高压透平入口压力
            data.P_lp_out - P_reheat,      % 再热压力必须大于低压透平出口压力
            P_intercool - data.P_hp_in,     % 中间冷却压力必须小于高压透平入口压力
            data.P_lp_out - P_intercool    % 中间冷却压力必须大于低压透平出口压力
        ];
    catch
        % 如果计算失败，给予惩罚值
        objectives = [100, 10]; % 最差的效率和成本
        constraints = [1, 1, 1, 1]; % 严重违反所有约束
    end
    
    % 恢复原始值
    assignin('base', 'P2', originalP2);
    assignin('base', 'P9', originalP9);
    assignin('base', 'alpha', originalAlpha);
    assignin('base', 'THTR', originalTHTR);
    assignin('base', 'TLTR', originalTLTR);
end

%% 目标函数1：负的热效率（最小化即最大化效率）
function f = objEfficiency(x, data)
    [~, objectives, ~] = evalSCO2Cycle(x, data);
    f = objectives(1);
end

%% 目标函数2：单位成本
function f = objCost(x, data)
    [~, objectives, ~] = evalSCO2Cycle(x, data);
    f = objectives(2);
end

%% 约束函数
function g = constraints(x, data)
    [~, ~, g] = evalSCO2Cycle(x, data);
end

%% 无效解修复函数
function x_repaired = repairSolution(x)
    % 确保所有变量在指定范围内
    lower = [10, 10, 0.2, 5, 5];
    upper = [20, 20, 0.5, 30, 30];
    x_repaired = max(lower, min(upper, x));
end

%% 调用PlatEMO求解
platemo('objFcn', {@objEfficiency, @objCost}, ...
        'encoding', [1, 1, 1, 1, 1], ... % 5个实数变量
        'lower', [10, 10, 0.2, 5, 5], ... % 变量下界
        'upper', [20, 20, 0.5, 30, 30], ... % 变量上界
        'conFcn', @constraints, ... % 约束函数
        'decFcn', @repairSolution, ... % 无效解修复函数
        'evalFcn', @evalSCO2Cycle, ... % 评价函数
        'data', data); % 问题数据
