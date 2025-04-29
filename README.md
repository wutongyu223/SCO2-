# 超临界CO2布雷顿循环多目标优化

此项目为超临界CO2布雷顿循环系统的多目标优化问题，使用PlatEMO平台上的NSGA-III算法进行优化。

## 问题描述

超临界CO2布雷顿循环是一种先进的能量转换循环，相比传统蒸汽循环具有更高的效率和更紧凑的尺寸。本项目定义了一个13点分流再压缩超临界CO2布雷顿循环模型，并对其关键参数进行多目标优化。

优化目标包括：
1. 最大化循环热效率
2. 最大化净输出功率
3. 最小化系统复杂度（由压力比和分流比组成的综合指标）

## 决策变量

本优化问题包含8个决策变量：
1. 高压透平入口压力(MPa) - 范围[18, 25]
2. 高压透平出口压力(MPa) - 范围[10, 18]
3. 低压透平出口压力(MPa) - 范围[5, 10]
4. 主压缩机a出口压力(MPa) - 范围[8, 15]
5. 分流比例 - 范围[0.2, 0.45]
6. 透平效率 - 范围[0.85, 0.95]
7. 主压缩机效率 - 范围[0.85, 0.95]
8. 回热器效率 - 范围[0.80, 0.95]

## 约束条件

优化过程中需满足以下约束条件：
1. 压力级别合适：P1 > P2 > P4
2. 主压缩机a出口压力小于高压透平入口压力：P8 < P1
3. 主压缩机a出口压力大于低压透平出口压力：P4 < P8

## 使用方法

### 前提条件
- 安装MATLAB (R2019b或更高版本)
- 安装PlatEMO平台
- 确保REFPROP功能可用（用于计算CO2热力学性质）

### 在PlatEMO中运行优化

1. 将`SCO2_BraytonCycle.m`文件复制到PlatEMO的`Problems/UserProblem`目录下。
2. 启动MATLAB，运行PlatEMO平台。
3. 在PlatEMO界面中，选择以下设置：
   - 算法：NSGA-III
   - 问题：SCO2_BraytonCycle
   - 种群大小：91（建议值，可调整）
   - 评价次数：10000（根据需要调整）
   - 目标数量：3

4. 点击"开始"按钮运行优化。

### 结果分析

优化完成后，PlatEMO将显示最优Pareto前沿，您可以：
1. 查看不同目标之间的权衡关系
2. 导出最优解的决策变量值和目标函数值
3. 通过调整PlatEMO的结果可视化选项，分析不同变量对目标的影响

## 代码结构

- `SCO2_BraytonCycle.m`：定义了PlatEMO可识别的多目标优化问题类
  - `CalObj`方法：计算目标函数值
  - `CalCon`方法：计算约束违反程度
  - `calcCyclePerformance`方法：执行循环热力学计算

## 注意事项

- 优化过程中可能出现无法收敛的情况，代码中通过try-catch结构进行了处理
- 为提高优化效率，计算过程简化为最多50次迭代
- 优化过程可能耗时较长，请耐心等待

## 参考文献

1. Dostal, V., Driscoll, M.J., Hejzlar, P., 2004. A supercritical carbon dioxide cycle for next generation nuclear reactors. Technical Report MIT-ANP-TR-100.
2. Lemmon, E.W., Huber, M.L., McLinden, M.O. NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties-REFPROP. 