%%------超临界CO2再压缩布雷顿循环优化程序------
% 主程序：使用PlatEMO平台的NSGA-III算法进行多目标优化

% 全局参考状态
REF.T0 = 298.15;   % 参考温度 [K]
REF.P0 = 101325;   % 参考压力 [Pa]
REF.h0 = refpropm('h','T',REF.T0,'P',REF.P0/1000,'CO2'); % [J/kg]
REF.s0 = refpropm('s','T',REF.T0,'P',REF.P0/1000,'CO2'); % [J/(kg·K)]

% 添加PlatEMO路径（根据实际安装路径修改）
% 获取脚本当前目录
scriptPath = fileparts(mfilename('fullpath'));
platemoPath = 'C:\Users\abc\Documents\GitHub\PlatEMO'; % 根据实际路径修改
addpath(genpath(platemoPath));

% 设置优化变量边界
% 变量定义: [P_reheat, P_intercool, x, THTR, TLTR]
lower = [1e7, 7e6, 0.2, 5, 5];     % 下界
upper = [3e7, 2e7, 0.5, 30, 30];   % 上界

% 定义目标函数（兼容PlatEMO接口 - 使用单元数组方式定义多个单目标函数）
f1 = @(x) efficiency_objective(x);   % 负热效率（最小化，等效于最大化效率）
f2 = @(x) cost_objective(x);         % 回热器成本

% 调用PlatEMO的NSGA-III算法进行优化
try    
    % 处理函数命名冲突 - 通过正确指定完整路径来执行内置函数
    % 而不是尝试移动文件
    disp('开始NSGA-III优化...');
    
    % 使用platemo函数调用NSGA-III算法
    [Dec, Obj] = platemo('algorithm', @NSGAIII, ...
                         'objFcn', {f1, f2}, ...  % 使用单元数组传入多个目标函数
                         'N', 100, ...            % 种群大小
                         'M', 2, ...              % 目标数量
                         'D', 5, ...              % 变量数量
                         'lower', lower, ...      % 变量下界
                         'upper', upper, ...      % 变量上界
                         'maxFE', 10000, ...      % 最大评价次数
                         'save', 10);             % 保存数量
    
    % 显示优化结果
    disp('优化完成！');
    disp(['最终种群大小: ', num2str(size(Dec, 1))]);
    
    % 绘制Pareto前沿
    figure;
    plot(-Obj(:,1), Obj(:,2), 'o');
    xlabel('热效率 η [%]');
    ylabel('单位成本 [USD/MW]');
    title('Pareto前沿 - NSGA-III优化结果(PlatEMO)');
    grid on;
    
    % 分析热效率最高的解
    [~, idx] = max(-Obj(:,1));
    best_solution = Dec(idx,:);
    disp('热效率最高的解:');
    disp(['P_reheat = ', num2str(best_solution(1)/1e6), ' MPa']);
    disp(['P_intercool = ', num2str(best_solution(2)/1e6), ' MPa']);
    disp(['分流比例 x = ', num2str(best_solution(3))]);
    disp(['THTR = ', num2str(best_solution(4)), ' K']);
    disp(['TLTR = ', num2str(best_solution(5)), ' K']);
    disp(['热效率 = ', num2str(-Obj(idx,1)), ' %']);
    disp(['单位成本 = ', num2str(Obj(idx,2)), ' 万美元/MW']);
    
    % 保存结果到文件
    result.Dec = Dec;
    result.Obj = Obj;
    result.best = best_solution;
    result.best_obj = Obj(idx,:);
    save('optimization_results_platemo.mat', 'result');
    
catch ME
    disp('优化过程中出错：');
    disp(ME.message);
    rethrow(ME);
end

%% 以下是PlatEMO优化所需的辅助函数

% 效率目标函数 - 返回负的热效率（需要最小化）
function f = efficiency_objective(x)
    try
        % 提取变量
        TLTR = x(5);
        P_reheat = x(1);
        P_intercool = x(2);
        alpha = x(3);
        THTR = x(4);
        
        % 检查输入范围
        if TLTR < 5 || TLTR > 30 || ...
           P_reheat < 1e7 || P_reheat > 3e7 || ...
           P_intercool < 7e6 || P_intercool > 2e7 || ...
           alpha < 0.2 || alpha > 0.5 || ...
           THTR < 5 || THTR > 30
            f = 1000; % 变量超出范围时返回惩罚值
            return;
        end
        
        % 调用热力学计算函数
        [eta_i, ~, ~] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
        
        % 返回负效率（最小化负效率等效于最大化效率）
        f = -eta_i * 100;  % 热效率 [%]
        
        % 确保返回的是有效实数
        if isnan(f) || isinf(f) || ~isreal(f) || f >= 0 || f < -100
            f = 1000;  % 无效效率
        end
    catch
        % 错误处理 - 返回一个默认的有效实数值
        f = 1000;
    end
    
    % 最后再次检查结果是否为有效实数
    if ~isfinite(f) || ~isreal(f)
        f = 1000;
    end
end

% 成本目标函数 - 返回回热器成本
function f = cost_objective(x)
    try
        % 提取变量
        TLTR = x(5);
        P_reheat = x(1);
        P_intercool = x(2);
        alpha = x(3);
        THTR = x(4);
        
        % 检查输入范围
        if TLTR < 5 || TLTR > 30 || ...
           P_reheat < 1e7 || P_reheat > 3e7 || ...
           P_intercool < 7e6 || P_intercool > 2e7 || ...
           alpha < 0.2 || alpha > 0.5 || ...
           THTR < 5 || THTR > 30
            f = 1e6; % 变量超出范围时返回惩罚值
            return;
        end
        
        % 调用热力学计算函数
        [~, CHTR_i, CLTR_i] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
        
        % 返回总成本
        f = (CHTR_i + CLTR_i) * 10^(-6);  % 总成本 [百万美元]
        
        % 确保返回的是有效实数
        if isnan(f) || isinf(f) || ~isreal(f) || f <= 0 || f > 1e6
            f = 1e6;  % 极高成本
        end
    catch
        % 错误处理 - 返回一个默认的有效实数值
        f = 1e6;
    end
    
    % 最后再次检查结果是否为有效实数
    if ~isfinite(f) || ~isreal(f)
        f = 1e6;
    end
end

% 热力学计算函数 - 基于原始cost2.m的核心计算逻辑
function [n, CHTR, CLTR] = calculate_cycle(TLTR, P_reheat, P_intercool, x, THTR)
% 预设默认返回值，以防计算失败
n = NaN;
CHTR = NaN;
CLTR = NaN;

try
%%---透平
n2 = 0.93;%透平效率
    P14 = P_reheat/1e6;   %进口压力，单位Mpa
T14 = str2double('840');  %进口温度，单位K
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
    P1 = P_intercool/1e6;  %出口压力，单位Mpa
S1 = S14;
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');
H1 = H14 + (H1 - H14)*n2;   %点1实际的焓
T1 = refpropm('t','H',H1,'P',P1*1000,'CO2');
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');

%%---主压缩机（1-x分流，效率N1）
N1 = 0.89;
T5 = str2double('305'); %令主压缩机进口温度32（效率曲线）
P5 = P1;                %主压缩机进口压力，P5等于透平出口压力P1
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
H5 = refpropm('h','T',T5,'P',P5*1000,'CO2');
S6 = S5;                %假设压缩完熵不变
P6 = P14;               %冷压缩机出口与透平进口相连，压力相同
H6 = refpropm('h','P',P6*1000,'S',S6,'CO2');
H6 = H5 + (H6 - H5)/N1;
T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');

%%---副压缩机（x分流,效率N2）
N2 = 0.89;
P9 = P1;                 %副压缩机进口与透平出口相连，压力相等

%%---冷却器（1-x分流）
P4 = P1;%冷却器进口与透平出口相连，压力相等

%%低温回热器
T7 = T6;
P7 = P6;
H7 = H6;
P8 = P7;

%%高温回热器
P10 = P14;
P11 = P10;
P3 = P1;
T2 = T1;
H2 = refpropm('h','T',T2,'P',P1*1000,'CO2');
P12 = P11;

%%---求解
e1 = 0.86;%低温回热器
e2 = 0.86;%高温回热器
n1 = 0.94;%锅炉效率
H111 = 300000;%初始假设点11的焓值

    object = 1;
    iter_count = 0;
    max_iter = 100; % 设置最大迭代次数，防止死循环

    while object && (iter_count < max_iter)
        iter_count = iter_count + 1;
        
        % 使用try-catch包装每次迭代中的计算
        try
        T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
        T3 = T11 + THTR;%高温回热器热端出口（利用端差计算）
        H3 = refpropm('H','T',T3,'P',P3*1000,'CO2');
        H12 = H2 - H3 + H111;
        T12 = refpropm('t','H',H12,'P',P12*1000,'CO2');
        T3 = refpropm('t','H',H3,'P',P3*1000,'CO2');
        T8 = T3 - TLTR;
        H8 = refpropm('h','T',T8,'P',P8*1000,'CO2');
        H9 = H3 - (1-x)*(H8 - H7);
        T9 = refpropm('t','H',H9,'P',P9*1000,'CO2');
        S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');
        S10 = S9;%        
        P10 = P14;
        H10 = refpropm('h','P',P10*1000,'S',S10,'CO2');
        H10 = H9 + (H10 - H9)/N2;%热压缩机实际出口焓值
        T10 = refpropm('t','P',P10*1000,'S',S10,'CO2');
        %冷却器
        H4 = H9;
        %分流汇合
        H11 = x*H10 + (1-x)*H8;
            
            % 检查H11/H111是否为有限值
            if ~isfinite(H11) || ~isfinite(H111) || H111 == 0
                return; % 如果无效，退出函数
            end
            
        if (0.98 <= H11/H111) && (H11/H111<= 1.02)
            object = 0;
        else
            H111 = H111 + 50;
        end
        catch
            % 如果计算过程中出错，直接返回默认值
            return;
        end
    end

    % 检查是否达到最大迭代次数
    if iter_count >= max_iter
        % 如果迭代未收敛，返回默认值
        return;
    end

    try
    n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1);%循环功：透平-压缩机耗功
        
        % 检查计算出的效率是否在合理范围
        if isnan(n) || ~isfinite(n) || n <= 0 || n > 1
            return; % 效率超出范围时返回默认值
        end
        
    W = 300;%输出净功率
    qm = (W*10^(6))/((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5));%(忽略发电机效率)-计算流量qm，单位kg/s
        
        % 检查流量是否为有限正值
        if isnan(qm) || ~isfinite(qm) || qm <= 0
            return; % 流量无效时返回默认值
        end
        
    W1 = x*(H10 - H9)*qm;
    W2 = (1-x)*(H6 - H5)*qm;
    %高温回热器
    e = exp(1);
    Qh = (H2 - H3)*qm;%单位W
    T212 = T2 - T12;
    T311 = T3 - T11;
        
        % 确保温差为正值，防止对数计算出错
        if T212 <= 0 || T311 <= 0 || max(T212,T311) <= 0 || min(T212,T311) <= 0
            return;
        end
        
    TH = (max(T212,T311) - min(T212,T311))/log(max(T212,T311)/min(T212,T311));
    UH = Qh/TH;
        
        % 检查UH是否为有限正值
        if isnan(UH) || ~isfinite(UH) || UH <= 0
            return;
        end
        
    ftHTR = 1;%最高温度小于550℃，温度修正系数为1
    fpHTR = 1;%压力修正系数
    aHTR = 49.45;%
    bHTR = 0.7544;%
    CHTR = aHTR*((UH)^bHTR)*ftHTR*fpHTR;
        
    %低温回热器
    T38 = T3 - T8;
    T47 = T9 - T7;
    Ql = (H3 - H9)*qm;%单位W
        
        % 确保温差为正值，防止对数计算出错
        if T38 <= 0 || T47 <= 0 || max(T38,T47) <= 0 || min(T38,T47) <= 0
            return;
        end
        
    TL = (max(T38,T47) - min(T38,T47))/log(max(T38,T47)/min(T38,T47));
        UL = Ql/TL;
        
        % 检查UL是否为有限正值
        if isnan(UL) || ~isfinite(UL) || UL <= 0
            return;
        end
        
    ftLTR = 1;%最高温度小于550℃，温度修正系数为1
    fpLTR = 1;%压力修正系数
    aLTR = 49.45;%
    bLTR = 0.7544;
    CLTR = aLTR*((UL)^bLTR)*ftLTR*fpLTR;

        % 最后检查计算结果是否为有限值
        if isnan(CHTR) || ~isfinite(CHTR) || CHTR <= 0 || ...
           isnan(CLTR) || ~isfinite(CLTR) || CLTR <= 0
            return;
        end
    catch
        % 如果计算过程中出错，返回默认值
        return;
    end
catch
    % 如果整个计算过程出错，返回默认值
    return;
end
end