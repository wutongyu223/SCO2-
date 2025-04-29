# while循环

在超临界CO2布雷顿循环计算中，`while`循环是一种重要的迭代求解工具，主要用于解决热平衡方程、收敛性计算和系统优化。

## 基本语法

```matlab
while 条件
    % 循环体内的代码
end
```

循环会一直执行直到条件表达式的结果为`false`或者遇到`break`语句。

## 在热力学计算中的应用

### 1. 求解热平衡方程

热平衡方程常常是隐式的，无法直接求解，需要通过迭代法逐步逼近解。

#### [[delete2.m代码结构]]中的应用
```matlab
object = 1;
while object
    T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
    T3 = T11 + THTR;%高温回热器热端出口（利用端差计算）
    H3 = refpropm('H','T',T3,'P',P3*1000,'CO2');
    ...
    H11 = x*H10 + (1-x)*H8;
    if (0.999 <= H11/H111) && (H11/H111<= 1.001)
        object = 0;
    else
        H111 = H111 + 50;
    end
end
```
- 从一个初始猜测值开始（`H111 = 300000`）
- 计算所有相关状态点
- 检查新计算的H11与假设值H111的比值是否接近1
- 如果满足收敛条件，退出循环；否则更新假设值继续迭代

### 2. 热平衡迭代优化

#### [[calautatlion.m代码结构]]中的应用
```matlab
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
    
    ...

    %检查收敛性
    dH5 = abs(H5 - H5_old)/abs(H5);
    dH6 = abs(H6 - H6_old)/abs(H6);
    dH15 = abs(H15 - H15_old)/abs(H15);
    dH16 = abs(H16 - H16_old)/abs(H16);
    
    if max([dH5, dH6, dH15, dH16]) < 0.001
        object = 0;
    end
end
```
- 引入迭代计数器`iter`和最大迭代次数`max_iter`防止无限循环
- 保存上一轮计算值用于计算变化率
- 使用松弛系数平滑收敛过程
- 基于多个变量的相对变化判断整体收敛性

### 3. 优化计算中的容错处理

#### [[cost2.m代码结构]]中的应用
```matlab
object = 1;
iter_count = 0;
max_iter = 100; % 设置最大迭代次数，防止死循环

while object && (iter_count < max_iter)
    iter_count = iter_count + 1;
    
    % 使用try-catch包装每次迭代中的计算
    try
        T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
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
- 结合[[try-catch结构]]处理迭代计算中可能出现的错误
- 放宽收敛判据（从0.999-1.001到0.98-1.02）增强优化过程的稳定性
- 设置最大迭代次数避免死循环

## 迭代收敛技巧

### 1. 松弛系数
为避免计算不稳定，常引入松弛系数α（通常0<α<1）来平滑变量更新：
```matlab
relax = 0.5; %松弛系数
H5 = H5_old + relax*(H5_new - H5_old);
```

### 2. 收敛判据
可以基于绝对误差或相对误差设计收敛判据：
```matlab
% 绝对误差
if abs(H11 - H111) < tolerance
    object = 0;
end

% 相对误差
if abs(H5 - H5_old)/abs(H5) < 0.001
    object = 0;
end
```

### 3. 最大迭代次数
防止无限循环，总是设置最大迭代次数：
```matlab
max_iter = 50;
iter = 0;
while condition && iter < max_iter
    iter = iter + 1;
    ...
end

if iter >= max_iter
    disp('警告：迭代未收敛，已达最大迭代次数');
end
```

## 常见问题及解决方案

### 1. 迭代发散
- 问题：迭代过程不收敛，变量值无限增大或震荡
- 解决方案：
  - 增加松弛系数
  - 选择更合理的初始值
  - 检查物理模型的合理性

### 2. 迭代收敛到物理上不合理的解
- 问题：计算收敛但结果不符合物理规律
- 解决方案：
  - 增加物理约束条件
  - 检查边界条件设置
  - 验证收敛解的物理合理性

### 3. 计算过程中出现数值问题
- 问题：出现NaN、Inf等无效值
- 解决方案：
  - 使用[[try-catch结构]]捕获错误
  - 检查函数输入参数范围
  - 增加数值有效性检查

## 与其他函数的关联
- [[refpropm函数]] - 在迭代中反复调用计算热力学性质
- [[max函数]] - 评估多变量收敛性
- [[abs函数]] - 计算误差大小
- [[try-catch结构]] - 处理迭代中的异常
- [[if语句]] - 设置收敛条件 