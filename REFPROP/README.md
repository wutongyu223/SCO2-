# 超临界CO₂布雷顿循环系统模拟与优化

本项目用于模拟和优化再热+分流再压缩(主压缩机中间冷却)+双透平结构的超临界CO₂布雷顿循环系统。项目使用MATLAB开发，包含热力学计算、㶲分析、经济分析和多目标优化。

## 系统结构

系统由以下主要组件构成：
- 双透平系统（高压透平和低压透平）
- 再热器
- 分流结构（部分流体通过副压缩机，部分通过主压缩机）
- 主压缩机（分为两级，中间有冷却）
- 副压缩机
- 高温回热器和低温回热器
- 冷却器和中间冷却器
- 加热器

## 系统状态点编号

1. 高压透平入口
2. 高压透平出口/再热器入口
3. 再热器出口/低压透平入口
4. 低压透平出口/高温回热器热侧入口
5. 高温回热器热侧出口/低温回热器热侧入口
6. 低温回热器热侧出口/分流点
7. 分流后主路：冷却器入口
8. 冷却器出口/主压缩机a入口
9. 主压缩机a出口/中间冷却入口
10. 中间冷却出口/主压缩机b入口
11. 主压缩机b出口/合流点
12. 分流后副路：副压缩机入口
13. 副压缩机出口/合流点
14. 合流点出口/低温回热器冷侧入口
15. 低温回热器冷侧出口/高温回热器冷侧入口
16. 高温回热器冷侧出口/加热器入口
17. 加热器出口/高压透平入口(同点1)

## 项目文件结构

项目包含以下MATLAB文件：

1. `main.m` - 主程序文件
   - 设定系统参数
   - 调用各个模块进行计算
   - 展示和保存结果

2. `calculate_cycle.m` - 循环计算模块
   - 计算所有状态点的热力学参数（温度、压力、焓、熵等）
   - 计算循环性能指标（热效率、净功率等）
   - 实现回热器迭代计算

3. `exergy_analysis.m` - 㶲分析模块
   - 计算各状态点的㶲
   - 计算各组件的㶲损失和㶲效率
   - 计算系统总㶲效率

4. `economic_analysis.m` - 经济分析模块
   - 计算各组件成本
   - 计算系统总投资成本
   - 计算单位功率成本（$/kW）和发电成本（$/kWh）

5. `optimize_cycle.m` - 优化模块
   - 使用NSGA-II算法进行多目标优化
   - 优化参数：再热压力、中间冷却压力、分流比例、回热器端差
   - 目标：最大热效率和最小系统成本

6. `plot_results.m` - 可视化模块
   - 绘制T-s图
   - 绘制优化参数与目标函数的关系图
   - 绘制帕累托前沿面

7. `sensitivity_analysis.m` - 敏感性分析模块
   - 分析参数变化对系统性能的影响
   - 生成敏感性分析图表

## 系统参数

系统主要参数及其典型值：

- 最高压力：约25 MPa
- 最低压力：约7.5 MPa
- 最高温度：约840 K (570℃)
- 最低温度：约305 K (32℃)
- 高压透平效率：约93%
- 低压透平效率：约93%
- 主压缩机效率：约89%
- 副压缩机效率：约89%
- 高温回热器效率：约86%
- 低温回热器效率：约86%
- 加热器效率：约94%

## 优化参数

需要优化的参数及其取值范围：

- 再热压力：约10-20 MPa
- 中间冷却压力：约10-20 MPa
- 分流比例：约0.2-0.5
- 高温回热器端差：约5-30 K
- 低温回热器端差：约5-30 K

## 技术要点

1. **热力学性质计算**
   - 使用REFPROP通过refpropm函数获取CO₂的热力学性质
   - 需要确保MATLAB能正确调用REFPROP

2. **回热器计算**
   - 实现迭代计算处理回热器两侧状态相互依赖问题
   - 设置合理的收敛准则和最大迭代次数

3. **能量平衡**
   - 确保系统满足能量守恒定律
   - 验证各组件的能量平衡

4. **㶲分析**
   - 计算各状态点的物理㶲
   - 计算各组件的㶲损失
   - 确保㶲平衡：总㶲输入 = 㶲输出 + 㶲损失

5. **多目标优化**
   - 使用NSGA-II算法实现多目标优化
   - 平衡系统效率和经济性

6. **结果可视化**
   - 使用MATLAB绘图功能展示结果
   - 包括T-s图、优化结果和敏感性分析

## 实施步骤

1. 开发循环计算模块 `calculate_cycle.m`
2. 开发㶲分析模块 `exergy_analysis.m`
3. 开发经济分析模块 `economic_analysis.m`
4. 实现优化算法 `optimize_cycle.m`
5. 开发可视化模块 `plot_results.m`
6. 实现敏感性分析 `sensitivity_analysis.m`
7. 整合所有模块到主程序 `main.m`
8. 测试系统并进行参数调整
9. 运行优化并分析结果

## 注意事项

- 确保已正确安装和配置REFPROP
- 迭代计算可能需要合理初值和收敛控制
- 优化过程计算量大，需要考虑计算效率
- 确保计算结果的物理合理性
- 仔细验证模型的正确性

## 预期输出

1. 最优系统参数配置
2. 各状态点的温度、压力、焓、熵等参数
3. 系统热效率和㶲效率
4. 各组件的能量传递、㶲损失和成本
5. T-s图和优化参数关系图
6. 敏感性分析结果
7. 最优系统的经济分析结果