# calautatlion.m 代码结构

`calautatlion.m`是超临界CO2再热布雷顿循环的设计计算程序，进行热力学计算并绘制T-S图。

## 代码功能
- 计算超临界CO2再热布雷顿循环的热力学参数
- 使用迭代法求解系统热平衡
- 计算各组件性能和循环热效率
- 绘制循环T-S图

## 代码结构

### 1. 高压透平参数设置 (1-13行)
```matlab
%%------初始参数
%---高压透平
n_hp_turb = 0.93; %高压透平效率
P1 = str2double('21'); %高压透平进口压力21，单位MPa
T1 = str2double('873'); %高压透平进口温度，873K;600℃
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');
```
- 使用[[str2double函数]]设置高压透平入口参数
- 使用[[refpropm函数]]计算热力学状态

### 2. 高压透平和再热器计算 (15-25行)
```matlab
%高压透平出口/再热器入口
P2 = str2double('15'); %高压透平出口压力，单位MPa
S2_is = S1; %等熵过程
H2_is = refpropm('h','P',P2*1000,'S',S2_is,'CO2'); %等熵焓变
H2 = H1 - (H1 - H2_is)*n_hp_turb; %实际焓变
T2 = refpropm('t','P',P2*1000,'H',H2,'CO2'); %点2的温度
S2 = refpropm('s','T',T2,'P',P2*1000,'CO2'); %点2实际的熵

%---再热器出口/低压透平入口
P3 = P2; %再热器压力不变
T3 = T1; %再热到与初始温度相同
H3 = refpropm('h','T',T3,'P',P3*1000,'CO2');
S3 = refpropm('s','T',T3,'P',P3*1000,'CO2');
```
- 基于等熵效率计算高压透平实际膨胀过程
- 设置再热器出口温度与初始温度相同

### 3. 低压透平计算 (27-35行)
```matlab
%---低压透平
n_lp_turb = 0.93; %低压透平效率
P4 = str2double('7.5729'); %低压透平出口压力，单位MPa
S4_is = S3; %等熵过程
H4_is = refpropm('h','P',P4*1000,'S',S4_is,'CO2'); %等熵焓变
H4 = H3 - (H3 - H4_is)*n_lp_turb; %实际焓变
T4 = refpropm('t','P',P4*1000,'H',H4,'CO2'); %点4的温度
S4 = refpropm('s','T',T4,'P',P4*1000,'CO2'); %点4实际的熵
```
- 基于等熵效率计算低压透平实际膨胀过程
- 计算低压透平出口工质状态

### 4. 端差设置和初始化 (37-48行)
```matlab
%---端差设置
THTR = 20; %高温回热器热端温差
TLTR = 20; %低温回热器热端温差

%---初始化其他状态点
%先初始化冷却器出口和主压缩机
P5 = P4; %高温回热器热侧入口压力等于低压透平出口压力
P6 = P5; %低温回热器热侧入口压力
P7 = P6; %冷却器入口压力（主路）
P8 = P7; %冷却器出口压力
T8 = str2double('305'); %冷却器出口温度32℃
H8 = refpropm('h','T',T8,'P',P8*1000,'CO2');
S8 = refpropm('s','T',T8,'P',P8*1000,'CO2');
```
- 设置高温和低温回热器端差
- 初始化回热器各侧压力
- 设置冷却器出口温度

### 5. 主压缩机计算 (50-71行)
```matlab
%---主压缩机a
N_mc_a = 0.89; %主压缩机a效率
P9 = str2double('12'); %主压缩机a出口压力/中间冷却入口压力
S9_is = S8; %等熵压缩
H9_is = refpropm('h','P',P9*1000,'S',S9_is,'CO2');
H9 = H8 + (H9_is - H8)/N_mc_a; %实际焓变
T9 = refpropm('t','P',P9*1000,'H',H9,'CO2');
S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');

%---中间冷却器
P10 = P9; %中间冷却出口压力
T10 = str2double('305'); %中间冷却出口温度32℃
H10 = refpropm('h','T',T10,'P',P10*1000,'CO2');
S10 = refpropm('s','T',T10,'P',P10*1000,'CO2');

%---主压缩机b
N_mc_b = 0.89; %主压缩机b效率
P11 = P1; %主压缩机b出口压力等于高压透平入口压力
S11_is = S10; %等熵压缩
H11_is = refpropm('h','P',P11*1000,'S',S11_is,'CO2');
H11 = H10 + (H11_is - H10)/N_mc_b; %实际焓变
T11 = refpropm('t','P',P11*1000,'H',H11,'CO2');
S11 = refpropm('s','T',T11,'P',P11*1000,'CO2');
```
- 两级压缩过程计算：主压缩机a和主压缩机b
- 中间设置中间冷却器
- 基于等熵效率计算实际压缩功耗

### 6. 分流设置和初始化 (73-93行)
```matlab
%---分流比例
alpha = 0.3333; %分流比例

%---合流点(14)
P14 = P11; %合流点压力等于主压缩机b出口压力
%由于还没确定H13，先暂时让H14 = H11
H14 = H11;
T14 = refpropm('t','P',P14*1000,'H',H14,'CO2');
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');

%---从合流点开始初步计算
P15 = P14; %低温回热器冷侧入口压力
P16 = P15; %高温回热器冷侧入口压力
P17 = P16; %加热器出口压力

%加热器出口温度等于高压透平入口温度
T17 = T1;
H17 = H1;
S17 = S1;
```
- 设置分流比例（0.3333）
- 初始化合流点参数
- 设置合流点后各状态点参数

### 7. 回热器初始计算 (95-117行)
```matlab
%---计算高温回热器冷侧出口
%初步估计高温回热器冷端温差
T16 = T4 - THTR;
H16 = refpropm('h','T',T16,'P',P16*1000,'CO2');
S16 = refpropm('s','T',T16,'P',P16*1000,'CO2');

%---初步计算低温回热器
T15 = T5 - TLTR; %低温回热器冷侧出口温度（根据端差）
H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');

%---计算高温回热器热侧出口
%我们需要基于能量平衡来计算高温回热器热侧出口
%此时可以使用H16和H15之间的能量差近似
Q_HT_cold = H16 - H15; %热量转移到冷侧
H5 = H4 - Q_HT_cold; %高温回热器热侧出口焓
T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');

%---计算低温回热器热侧出口
Q_LT_cold = H15 - H14; %热量转移到冷侧
H6 = H5 - Q_LT_cold; %低温回热器热侧出口焓
T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
```
- 基于设定的端差计算高温和低温回热器冷侧出口温度
- 基于能量平衡初步估计回热器热侧出口状态

### 8. 副压缩机计算 (119-136行)
```matlab
%---分流点
%分流后每一路的热力学状态相同
P12 = P6; %副压缩机入口压力（副路）
T7 = T6;
T12 = T6;
H7 = H6;
H12 = H6;
S7 = S6;
S12 = S6;

%---副压缩机
N_rc = 0.89; %副压缩机效率
P13 = P1; %副压缩机出口压力等于高压透平入口压力
S13_is = S12; %等熵压缩
H13_is = refpropm('h','P',P13*1000,'S',S13_is,'CO2');
H13 = H12 + (H13_is - H12)/N_rc; %实际焓变
T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
```
- 设置分流点后主路和副路的状态相同
- 计算副压缩机的压缩过程
- 基于等熵效率计算实际功耗

### 9. 合流点与热平衡迭代 (138-224行)
```matlab
%---更新合流点计算
H14 = (1-alpha)*H11 + alpha*H13; %能量守恒，质量加权平均
T14 = refpropm('t','P',P14*1000,'H',H14,'CO2');
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');

%%---热平衡迭代求解
e1 = 0.86; %低温回热器效率
e2 = 0.86; %高温回热器效率
n_heater = 0.94; %加热器效率

%热平衡迭代求解
object = 1;
iter = 0;
max_iter = 50;

while object && iter < max_iter
    iter = iter + 1;
    
    %保存上一轮值用于比较收敛
    H5_old = H5;
    H6_old = H6;
    H15_old = H15;
    H16_old = H16;
    
    %更新高温回热器
    Q_HT_hot = H4 - H5; %高温回热器热侧放热
    Q_HT_cold = H16 - H15; %高温回热器冷侧吸热
    
    %检查能量守恒并进行调整
    dQ_HT = Q_HT_hot - Q_HT_cold;
    %调整高温回热器热侧出口焓
    H5_new = H4 - Q_HT_cold; %高温回热器热侧出口焓
    
    %带松弛系数的更新
    relax = 0.5; %松弛系数
    H5 = H5_old + relax*(H5_new - H5_old);
    T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
    S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
    
    ...

    %检查收敛性
    dH5 = abs(H5 - H5_old)/abs(H5);
    dH6 = abs(H6 - H6_old)/abs(H6);
    dH15 = abs(H15 - H15_old)/abs(H15);
    dH16 = abs(H16 - H16_old)/abs(H16);
    
    if max([dH5, dH6, dH15, dH16]) < 0.001
        object = 0;
    end
```
- 计算合流点参数：质量加权平均焓值
- 使用[[while循环]]迭代求解系统热平衡
- 添加松弛系数提高收敛稳定性
- 使用相对误差判断收敛条件

### 10. 系统性能计算 (226-250行)
```matlab
%---计算循环效率
%功率计算
W_hp_turb = H1 - H2; %高压透平做功
W_lp_turb = H3 - H4; %低压透平做功
W_mc_a = H9 - H8; %主压缩机a耗功
W_mc_b = H11 - H10; %主压缩机b耗功
W_rc = H13 - H12; %副压缩机耗功
W_comp_total = (1-alpha)*(W_mc_a + W_mc_b) + alpha*W_rc; %总压缩功率
W_turb_total = W_hp_turb + W_lp_turb; %总透平功率
W_net = W_turb_total - W_comp_total; %净输出功率

%热量计算
Q_heater = H17 - H16; %加热器热量
Q_reheater = H3 - H2; %再热器热量
Q_cooler = (1-alpha)*(H7 - H8); %主路冷却器热量
Q_intercooler = (1-alpha)*(H9 - H10); %中间冷却器热量
Q_in_total = Q_heater + Q_reheater; %总输入热量

%循环热效率
eta = W_net / Q_in_total;
```
- 分别计算高低压透平做功、各压缩机耗功
- 计算总功率、净功率、热量和循环热效率
- 考虑分流比例的影响

### 11. 结果输出与T-S图绘制 (252-402行)
```matlab
%输出结果
disp(['高压透平做功: ' num2str(W_hp_turb) ' J/kg'])
disp(['低压透平做功: ' num2str(W_lp_turb) ' J/kg'])
disp(['总压缩功率: ' num2str(W_comp_total) ' J/kg'])
disp(['净输出功率: ' num2str(W_net) ' J/kg'])
disp(['总输入热量: ' num2str(Q_in_total) ' J/kg'])
disp(['循环热效率: ' num2str(eta*100) ' %'])

% 调用T-S绘图函数
plot_ts_diagram(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17,...
                S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17,...
                P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17,...
                alpha, W_net, eta);
```
- 使用[[disp函数]]和[[num2str函数]]输出计算结果
- 调用[[plot_ts_diagram函数]]绘制T-S图

### 12. T-S图绘制函数 (261-339行)
```matlab
function plot_ts_diagram(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17,...
                         S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15, S16, S17,...
                         P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17,...
                         alpha, W_net, eta)
    % 创建一个新图
    figure;
    hold on;
    
    % 设置标题和标签
    title('超临界CO₂布雷顿循环T-S图');
    xlabel('熵 s [J/(kg·K)]');
    ylabel('温度 T [K]');
    
    % 主循环连线
    % 高压透平 (1→2)
    line1 = plot([S1,S2], [T1,T2], 'b-', 'LineWidth', 1.5);
    
    % 再热器 (2→3) - 等压加热
    line2 = plot([S2,S3], [T2,T3], 'r-', 'LineWidth', 1.5);
    
    ...
    
    % 标记主循环状态点
    main_points = [1,2,3,4,5,6,7,8,9,10,11,14,15,16,17];
    main_T = [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T14,T15,T16,T17];
    main_S = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S14,S15,S16,S17];
    
    for i = 1:length(main_points)
        plot(main_S(i), main_T(i), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        text(main_S(i), main_T(i), [' ', num2str(main_points(i))], 'FontSize', 8);
    end
    
    ...
    
    % 添加图例
    legend([line1, line2, line16, plot(NaN,NaN,'ro','MarkerFaceColor','r'), plot(NaN,NaN,'bo','MarkerFaceColor','b')], ...
           {'主循环路径', '吸热过程', '分流路径', '主循环状态点', '分流状态点'}, ...
           'Location', 'Best', 'FontSize', 8);
```
- 使用[[figure函数]]创建新图
- 使用[[hold函数]]保持图形状态
- 使用[[plot函数]]绘制各状态点和连线
- 使用[[text函数]]标注状态点
- 使用[[legend函数]]添加图例

## 关键参数

| 参数 | 值 | 描述 |
|------|-----|------|
| n_hp_turb | 0.93 | 高压透平效率 |
| n_lp_turb | 0.93 | 低压透平效率 |
| N_mc_a | 0.89 | 主压缩机a效率 |
| N_mc_b | 0.89 | 主压缩机b效率 |
| N_rc | 0.89 | 副压缩机效率 |
| e1 | 0.86 | 低温回热器效率 |
| e2 | 0.86 | 高温回热器效率 |
| THTR | 20 | 高温回热器端差 |
| TLTR | 20 | 低温回热器端差 |
| alpha | 0.3333 | 分流比 |
| relax | 0.5 | 松弛系数 |

## 程序特点
1. 使用迭代法求解热平衡，确保能量守恒
2. 采用松弛系数提高收敛稳定性
3. 详细计算各组件性能指标
4. 绘制完整T-S图，直观展示循环过程
5. 考虑再热和分流对循环效率的影响

## 相关函数
- [[refpropm函数]] - 计算CO2热力学性质
- [[str2double函数]] - 参数初始化
- [[while循环]] - 迭代求解热平衡
- [[plot函数]] - 绘制T-S图
- [[disp函数]] - 显示结果
- [[num2str函数]] - 数值转字符串
- [[figure函数]] - 创建图形窗口
- [[hold函数]] - 保持图形状态
- [[abs函数]] - 计算误差绝对值
- [[max函数]] - 判断收敛性

## 相关代码
- [[delete2.m代码结构]] - 基础再压缩循环计算
- [[cost2.m代码结构]] - 循环优化计算 