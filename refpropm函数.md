# refpropm函数

`refpropm`是REFPROP（REFerence fluid PROPerties）库的MATLAB接口函数，用于计算各种工质在不同状态下的热力学性质。在超临界CO2布雷顿循环计算中，它是热力学状态计算的核心函数。

## 函数语法

```matlab
result = refpropm(output_property, input_property1, input_value1, input_property2, input_value2, fluid)
```

## 参数说明

### 输出性质代码 (output_property)
- `'T'` - 温度 [K]
- `'P'` - 压力 [kPa]
- `'D'` - 密度 [kg/m³]
- `'H'` - 比焓 [J/kg]
- `'S'` - 比熵 [J/(kg·K)]
- `'U'` - 内能 [J/kg]
- `'C'` - 定压比热容 [J/(kg·K)]
- `'O'` - 定容比热容 [J/(kg·K)]
- `'V'` - 声速 [m/s]

### 输入性质参数 (input_property)
与输出性质代码相同，用任意两个已知参数确定热力学状态

### 输入数值 (input_value)
与输入性质参数对应的数值

### 工质名称 (fluid)
可以是单一工质或混合工质，例如：`'CO2'`, `'water'`, `'R134a'`等

## 在代码中的应用示例

### 1. 基于温度和压力计算熵值
```matlab
% 计算CO2在873K，20MPa下的熵值
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
```

### 2. 基于压力和熵值计算焓值（等熵过程）
```matlab
% 计算CO2在等熵膨胀到7.6MPa过程中的焓值
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');
```

### 3. 基于压力和焓值计算温度
```matlab
% 计算CO2在7.6MPa，已知焓值H1下的温度
T1 = refpropm('t','P',P1*1000,'H',H1,'CO2');
```

## 注意事项

### 压力单位转换
`refpropm`函数要求压力单位为kPa，而代码中常用MPa，因此需要进行单位转换：
```matlab
P_kPa = P_MPa * 1000;
```

### 常见计算组合
1. 等熵过程计算：已知起始点温度和压力，以及终点压力
   ```matlab
   S_start = refpropm('s','T',T_start,'P',P_start*1000,'CO2');
   H_end_isentropic = refpropm('h','P',P_end*1000,'S',S_start,'CO2');
   ```

2. 等压过程计算：已知压力和两个温度
   ```matlab
   H_start = refpropm('h','T',T_start,'P',P*1000,'CO2');
   H_end = refpropm('h','T',T_end,'P',P*1000,'CO2');
   Q = H_end - H_start; % 焓变 = 吸热量
   ```

3. 实际膨胀/压缩过程计算：使用等熵效率修正
   ```matlab
   % 膨胀过程（透平）
   H_actual = H_start - (H_start - H_end_isentropic) * eta_turbine;
   
   % 压缩过程（压缩机）
   H_actual = H_start + (H_end_isentropic - H_start) / eta_compressor;
   ```

## 在各文件中的应用

### [[delete2.m代码结构]]
- 计算透平膨胀过程
- 计算压缩机压缩过程
- 迭代计算回热器热平衡

### [[cost2.m代码结构]]
- 优化变量条件下的热力学计算
- 在[[calculate_cycle函数]]中广泛使用

### [[calautatlion.m代码结构]]
- 计算再热布雷顿循环的各状态点
- 迭代求解热平衡方程

## 常见错误处理

当传入的参数组合不存在物理解或超出工作范围时，函数可能返回NaN或抛出错误。在代码中通常采用以下方式处理：

```matlab
try
    H = refpropm('h','T',T,'P',P*1000,'CO2');
    if isnan(H) || ~isfinite(H)
        % 处理无效结果
    end
catch
    % 处理计算错误
end
```

## 与其他函数的关联
- [[str2double函数]] - 为refpropm提供输入参数
- [[while循环]] - 迭代使用refpropm求解热平衡
- [[try-catch结构]] - 处理refpropm计算错误
- [[calculate_cycle函数]] - 包含多次调用refpropm

## 相关资源
- REFPROP官方文档：https://www.nist.gov/srd/refprop
- 需要REFPROP软件安装并配置MATLAB接口 