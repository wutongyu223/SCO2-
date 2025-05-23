# delete2.m 代码结构

`delete2.m`是超临界CO2再压缩布雷顿循环基本热力计算程序，计算固定参数下循环的热力性能和能源品位。

## 代码功能
- 计算固定[[分流比]]下超临界CO2再压缩布雷顿循环的热力参数
- 计算循环热效率和系统能源品位
- 分析各组件的能量转换过程

## 代码结构

### 1. 初始参数设置 (1-36行)
```matlab
n2 = 0.93;%透平效率
P14 = str2double('20');   %进口压力20，单位Mpa
T14 = str2double('873');  %进口温度，873K;600℃
```
- 使用[[str2double函数]]设置基本参数
- 使用[[refpropm函数]]计算初始热力学状态
- 参数包括：透平效率、压缩机效率、回热器效率、分流比等

### 2. 透平计算 (1-12行)
```matlab
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
P1 = str2double('7.6');   %出口压力，单位Mpa
S1 = S14;
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');
T1 = refpropm('t','H',H1,'P',P1*1000,'CO2');%点1的温度
H1 = H14 + (H1 - H14)*n2; %点1实际的焓
S1 = refpropm('s','P',P1*1000,'H',H1,'CO2');%点1实际的熵
```
- 基于等熵效率计算实际膨胀过程
- 计算透平出口工质状态点

### 3. 压缩机计算 (14-37行)
```matlab
N1 = 0.89;
T5 = str2double('305'); %令主压缩机进口温度32℃
P5 = P1;                %主压缩机进口压力，P5等于透平出口压力P1
...
H6 = H5 + (H6 - H5)/N1;%点6实际焓值
```
- 主压缩机（1-x分流）计算
- 副压缩机（x分流）计算
- 基于等熵效率计算实际压缩过程

### 4. 回热器初始设置 (38-54行)
```matlab
P7 = P6;
H7 = H6;
P8 = P7;
...
P10 = P14;
P11 = P10;
P3 = P1;
```
- 设置高温回热器与低温回热器初始参数
- 设置回热器端差和效率

### 5. 热平衡迭代求解 (55-93行)
```matlab
e1 = 0.86;%低温回热器
e2 = 0.86;%高温回热器
n1 = 0.94;%锅炉效率
H111 = 300000;%初始假设点11的焓值

THTR = 15;%高温回热器冷端端差
TLTR = 15;%低温回热器冷端端差

object = 1;
while object
    T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
    T3 = T11 + THTR;%高温回热器热端出口（利用端差计算）
    ...
    if (0.999 <= H11/H111) && (H11/H111<= 1.001)
        object = 0;
    else
        H111 = H111 + 50;
    end
end
```
- 使用[[while循环]]迭代求解系统热平衡
- 基于设定的端差计算回热器出入口状态
- 计算分流汇合点能量平衡
- 使用收敛条件判断迭代停止

### 6. 系统性能计算 (94-95行)
```matlab
n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1);
W1 = x*(H10 - H9) + (1-x)*(H6 - H5);
```
- 计算循环总热效率
- 计算压缩功耗

### 7. 能源品位分析 (96-148行)
```matlab
T0 = 25 + 273;
H0 = refpropm('h','T',T0,'P',101.325,'CO2');
S0 = refpropm('s','T',T0,'P',101.325,'CO2');
%---透平
E14 = (H14/1000 - H0/1000) - T0*(S14/1000 - S0/1000);
E1 = (H1/1000 - H0/1000) - T0*(S1/1000 - S0/1000);
ET = (E14 - E1) - (H14/1000 - H1/1000);
NT = 1-(ET)/(E14 - E1);
```
- 基于环境温度(25℃)计算各组件能源品位
- 分析透平、压缩机、回热器、冷却器的热力学性能
- 计算各组件的能量转换效率

### 8. 结果输出 (150-156行)
```matlab
disp(['分流比',num2str(x)])
disp(['THTR',num2str(THTR)])
disp(['TLTR',num2str(TLTR)])
disp(['循环效率',num2str(n)])
```
- 使用[[disp函数]]和[[num2str函数]]输出计算结果
- 包括分流比、端差、循环效率、能源品位等关键参数

## 关键参数

| 参数 | 值 | 描述 |
|------|-----|------|
| n2 | 0.93 | 透平效率 |
| N1 | 0.89 | 主压缩机效率 |
| N2 | 0.89 | 副压缩机效率 |
| e1 | 0.86 | 低温回热器效率 |
| e2 | 0.86 | 高温回热器效率 |
| n1 | 0.94 | 锅炉效率 |
| x | 0.35 | 分流比 |
| THTR | 15 | 高温回热器端差 |
| TLTR | 15 | 低温回热器端差 |

## 相关函数
- [[refpropm函数]] - 计算CO2热力学状态
- [[str2double函数]] - 参数初始化
- [[disp函数]] - 显示结果
- [[num2str函数]] - 数值转字符串

## 相关代码
- [[cost2.m代码结构]] - 基于本代码的优化版本
- [[calautatlion.m代码结构]] - 类似循环的再热版本 