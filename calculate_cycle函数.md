# calculate_cycle函数

`calculate_cycle`是[[cost2.m代码结构]]中的核心热力计算函数，用于在给定优化变量下计算超临界CO2再压缩布雷顿循环的热效率和回热器成本。

## 函数定义

```matlab
function [n, CHTR, CLTR] = calculate_cycle(TLTR, P_reheat, P_intercool, x, THTR)
```

## 输入参数

| 参数 | 描述 | 单位 |
|------|------|------|
| TLTR | 低温回热器端差 | K |
| P_reheat | 再热压力 | Pa |
| P_intercool | 中间冷却压力 | Pa |
| x | 分流比 | - |
| THTR | 高温回热器端差 | K |

## 输出参数

| 参数 | 描述 | 单位 |
|------|------|------|
| n | 循环热效率 | - |
| CHTR | 高温回热器成本 | USD |
| CLTR | 低温回热器成本 | USD |

## 函数结构

### 1. 错误处理初始化
```matlab
% 预设默认返回值，以防计算失败
n = NaN;
CHTR = NaN;
CLTR = NaN;

try
```
- 设置默认返回值为NaN
- 使用[[try-catch结构]]包装整个函数体

### 2. 透平计算
```matlab
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
```
- 将输入的Pa单位压力转换为MPa
- 使用[[refpropm函数]]计算透平进出口状态
- 基于等熵效率计算实际膨胀过程

### 3. 压缩机计算
```matlab
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
```
- 主压缩机计算：根据等熵效率计算实际压缩过程
- 副压缩机使用相同的计算方法
- 冷却器参数设置

### 4. 热平衡迭代计算
```matlab
object = 1;
iter_count = 0;
max_iter = 100; % 设置最大迭代次数，防止死循环

while object && (iter_count < max_iter)
    iter_count = iter_count + 1;
    
    % 使用try-catch包装每次迭代中的计算
    try
    T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
    T3 = T11 + THTR;%高温回热器热端出口（利用端差计算）
    ...
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
```
- 使用[[while循环]]迭代求解系统热平衡
- 设置最大迭代次数避免死循环
- 迭代计算过程中嵌套[[try-catch结构]]处理计算错误
- 收敛判据为H11/H111比值在0.98到1.02之间

### 5. 性能计算
```matlab
n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1);%循环功：透平-压缩机耗功
    
    % 检查计算出的效率是否在合理范围
    if isnan(n) || ~isfinite(n) || n <= 0 || n > 1
        return; % 效率超出范围时返回默认值
    end
    
W = 300;%输出净功率
qm = (W*10^(6))/((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5));%(忽略发电机效率)-计算流量qm，单位kg/s
```
- 计算循环热效率
- 检查计算结果是否在物理合理范围内
- 基于设定净功率计算工质流量

### 6. 回热器成本计算
```matlab
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
```
- 计算高温回热器的对数平均温差
- 计算传热系数-面积乘积（UA值）
- 使用成本估算公式计算回热器成本
- 低温回热器使用相同方法计算

### 7. 结果验证
```matlab
% 最后检查计算结果是否为有限值
if isnan(CHTR) || ~isfinite(CHTR) || CHTR <= 0 || ...
   isnan(CLTR) || ~isfinite(CLTR) || CLTR <= 0
    return;
end
```
- 检查计算得到的回热器成本是否为有效值
- 使用[[isnan函数]]和[[isfinite函数]]进行检查

## 函数的健壮性设计

`calculate_cycle`函数具有多层错误处理机制，确保在优化过程中不会因单点计算失败而导致整个优化中断：

1. **预设默认返回值**：函数开始设置默认返回值为NaN
2. **整体try-catch包装**：捕获所有可能的错误
3. **迭代计算的嵌套try-catch**：单次迭代失败时不影响整体计算
4. **多处物理合理性检查**：
   - 效率值在(0,1)范围内
   - 温差为正值
   - 传热系数为有限正值
   - 成本为有限正值
5. **最大迭代次数限制**：防止死循环
6. **宽松的收敛判据**：从0.999-1.001放宽到0.98-1.02

## 与优化目标函数的关系

`calculate_cycle`函数为[[cost2.m代码结构]]中的两个目标函数提供核心计算支持：

1. **效率目标函数**：
```matlab
function f = efficiency_objective(x)
    ...
    [eta_i, ~, ~] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
    f = -eta_i * 100;  % 热效率 [%]，取负值用于最小化
```

2. **成本目标函数**：
```matlab
function f = cost_objective(x)
    ...
    [~, CHTR_i, CLTR_i] = calculate_cycle(TLTR, P_reheat, P_intercool, alpha, THTR);
    f = (CHTR_i + CLTR_i) * 10^(-6);  % 总成本 [百万美元]
```

## 与其他代码的关系

`calculate_cycle`函数基本逻辑与[[delete2.m代码结构]]相似，但增加了以下改进：

1. 接受优化变量作为输入参数
2. 增强的错误处理机制
3. 增加了回热器成本计算
4. 更多的物理合理性检查
5. 返回优化所需的性能参数

## 相关函数
- [[refpropm函数]] - 热力学性质计算
- [[while循环]] - 迭代热平衡求解
- [[try-catch结构]] - 错误处理
- [[isnan函数]]和[[isfinite函数]] - 数值有效性检查
- [[log函数]] - 计算对数平均温差 