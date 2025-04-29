# cost2.m 代码结构

`cost2.m`是超临界CO2再压缩布雷顿循环的多目标优化程序，使用PlatEMO平台的NSGA-III算法寻找最优循环参数组合。

## 代码功能
- 优化超临界CO2再压缩布雷顿循环的设计参数
- 多目标优化：最大化热效率，最小化回热器成本
- 得到Pareto最优解集，分析不同权衡方案

## 代码结构

### 1. 全局参考状态设置 (3-7行)
```matlab
% 全局参考状态
REF.T0 = 298.15;   % 参考温度 [K]
REF.P0 = 101325;   % 参考压力 [Pa]
REF.h0 = refpropm('h','T',REF.T0,'P',REF.P0/1000,'CO2'); % [J/kg]
REF.s0 = refpropm('s','T',REF.T0,'P',REF.P0/1000,'CO2'); % [J/(kg·K)]
```
- 定义环境参考状态（温度、压力）
- 使用[[refpropm函数]]计算参考状态下的熵和焓

### 2. 添加优化平台路径 (9-13行)
```matlab
% 添加PlatEMO路径（根据实际安装路径修改）
% 获取脚本当前目录
scriptPath = fileparts(mfilename('fullpath'));
platemoPath = 'C:\Users\abc\Documents\GitHub\PlatEMO'; % 根据实际路径修改
addpath(genpath(platemoPath));
```
- 使用[[fileparts函数]]和[[mfilename函数]]获取当前脚本路径
- 使用[[addpath函数]]和[[genpath函数]]添加PlatEMO平台路径

### 3. 设置优化变量边界 (15-18行)
```matlab
% 变量定义: [P_reheat, P_intercool, x, THTR, TLTR]
lower = [1e7, 7e6, 0.2, 5, 5];     % 下界
upper = [3e7, 2e7, 0.5, 30, 30];   % 上界
```
- 设置5个优化变量的上下边界
- 分别为：再热压力、中间冷却压力、分流比、高温回热器端差、低温回热器端差

### 4. 定义目标函数 (20-23行)
```matlab
% 定义目标函数（兼容PlatEMO接口 - 使用单元数组方式定义多个单目标函数）
f1 = @(x) efficiency_objective(x);   % 负热效率（最小化，等效于最大化效率）
f2 = @(x) cost_objective(x);         % 回热器成本
```
- 定义两个目标函数的句柄
- 热效率目标使用负值便于最小化处理

### 5. NSGA-III算法调用 (25-51行)
```matlab
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
```
- 使用[[try-catch结构]]防止优化过程出错
- 通过[[platemo函数]]调用NSGA-III多目标优化算法
- 设置种群大小、目标数量、优化变量数量和最大评价次数

### 6. 优化结果可视化与保存 (53-79行)
```matlab
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
    ...
    
    % 保存结果到文件
    result.Dec = Dec;
    result.Obj = Obj;
    result.best = best_solution;
    result.best_obj = Obj(idx,:);
    save('optimization_results_platemo.mat', 'result');
```
- 使用[[figure函数]]和[[plot函数]]绘制Pareto前沿
- 使用[[max函数]]找出热效率最高的解
- 使用[[disp函数]]和[[num2str函数]]显示优化结果
- 使用[[save函数]]保存优化结果到文件

### 7. 效率目标函数 (84-120行)
```matlab
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
           ...
            f = 1000; % 变量超出范围时返回惩罚值
            return;
        end
        
        % 调用热力学计算函数
        [eta_i, ~, ~] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
        
        % 返回负效率（最小化负效率等效于最大化效率）
        f = -eta_i * 100;  % 热效率 [%]
```
- 定义效率目标函数，接收优化变量向量作为输入
- 使用[[try-catch结构]]处理计算错误
- 检查变量是否在有效范围内
- 调用[[calculate_cycle函数]]计算热效率
- 返回负的热效率作为优化目标（最小化）

### 8. 成本目标函数 (122-152行)
```matlab
function f = cost_objective(x)
    try
        % 提取变量
        TLTR = x(5);
        P_reheat = x(1);
        P_intercool = x(2);
        alpha = x(3);
        THTR = x(4);
        
        % 检查输入范围
        ...
        
        % 调用热力学计算函数
        [~, CHTR_i, CLTR_i] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
        
        % 返回总成本
        f = (CHTR_i + CLTR_i) * 10^(-6);  % 总成本 [百万美元]
```
- 定义成本目标函数，接收优化变量向量作为输入
- 使用[[try-catch结构]]处理计算错误
- 检查变量是否在有效范围内
- 调用[[calculate_cycle函数]]计算高温和低温回热器成本
- 返回总成本作为优化目标（最小化）

### 9. 热力计算函数 (154-356行)
```matlab
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
...
```
- 定义核心热力计算函数[[calculate_cycle函数]]
- 接收优化变量作为输入
- 计算布雷顿循环的热效率和回热器成本
- 包含多重错误检查和[[try-catch结构]]
- 基本计算逻辑与[[delete2.m代码结构]]类似，但针对优化环境进行了健壮性增强

## 优化变量

| 变量名 | 含义 | 下限 | 上限 | 单位 |
|--------|------|------|------|------|
| P_reheat | 再热压力 | 10 | 30 | MPa |
| P_intercool | 中间冷却压力 | 7 | 20 | MPa |
| alpha | 分流比 | 0.2 | 0.5 | - |
| THTR | 高温回热器端差 | 5 | 30 | K |
| TLTR | 低温回热器端差 | 5 | 30 | K |

## 目标函数

| 目标 | 函数名 | 优化方向 | 描述 |
|------|--------|----------|------|
| 热效率 | efficiency_objective | 最大化 | 系统热效率（代码中用负值表示） |
| 成本 | cost_objective | 最小化 | 高温和低温回热器的总成本 |

## 使用的PlatEMO参数

| 参数 | 值 | 描述 |
|------|-----|------|
| algorithm | NSGAIII | 多目标优化算法 |
| N | 100 | 种群大小 |
| M | 2 | 目标数量 |
| D | 5 | 优化变量数量 |
| maxFE | 10000 | 最大评价次数 |

## 相关函数
- [[platemo函数]] - 多目标优化算法调用
- [[refpropm函数]] - 热力学计算
- [[calculate_cycle函数]] - 核心热力计算
- [[plot函数]] - 绘制Pareto前沿
- [[save函数]] - 保存优化结果
- [[try-catch结构]] - 错误处理

## 相关代码
- [[delete2.m代码结构]] - 单一参数下的循环计算
- [[calautatlion.m代码结构]] - 再热循环设计计算 