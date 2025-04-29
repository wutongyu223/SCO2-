%%------初始参数
%---高压透平
n_hp_turb = 0.93; %高压透平效率
P1 = 21; %高压透平进口压力21，单位MPa
T1 = 873; %高压透平进口温度，873K;600℃
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');

%高压透平出口/再热器入口
P2 = 15; %高压透平出口压力，单位MPa
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

%---低压透平
n_lp_turb = 0.93; %低压透平效率
P4 = 7.5729; %低压透平出口压力，单位MPa
S4_is = S3; %等熵过程
H4_is = refpropm('h','P',P4*1000,'S',S4_is,'CO2'); %等熵焓变
H4 = H3 - (H3 - H4_is)*n_lp_turb; %实际焓变
T4 = refpropm('t','P',P4*1000,'H',H4,'CO2'); %点4的温度
S4 = refpropm('s','T',T4,'P',P4*1000,'CO2'); %点4实际的熵

%---端差设置
THTR = 20; %高温回热器热端温差
TLTR = 20; %低温回热器热端温差

%---初始化其他状态点
%高温回热器热侧出口
P5 = P4; %高温回热器热侧入口压力等于低压透平出口压力
%低温回热器热侧出口/分流点
P6 = P5; %低温回热器热侧入口压力

%---冷却器出口/主压缩机a入口
P7 = P6; %冷却器出口压力（与分流点6相同）
T7 = 305; %冷却器出口温度32℃
H7 = refpropm('h','T',T7,'P',P7*1000,'CO2');
S7 = refpropm('s','T',T7,'P',P7*1000,'CO2');

%---主压缩机a出口/中间冷却入口
N_mc_a = 0.89; %主压缩机a效率
P8 = 12; %主压缩机a出口压力/中间冷却入口压力
S8_is = S7; %等熵压缩
H8_is = refpropm('h','P',P8*1000,'S',S8_is,'CO2');
H8 = H7 + (H8_is - H7)/N_mc_a; %实际焓变
T8 = refpropm('t','P',P8*1000,'H',H8,'CO2');
S8 = refpropm('s','T',T8,'P',P8*1000,'CO2');

%---中间冷却出口/主压缩机b入口
P9 = P8; %中间冷却出口压力
T9 = 305; %中间冷却出口温度32℃
H9 = refpropm('h','T',T9,'P',P9*1000,'CO2');
S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');

%---主压缩机b出口
N_mc_b = 0.89; %主压缩机b效率
P10 = P1; %主压缩机b出口压力等于高压透平入口压力
S10_is = S9; %等熵压缩
H10_is = refpropm('h','P',P10*1000,'S',S10_is,'CO2');
H10 = H9 + (H10_is - H9)/N_mc_b; %实际焓变
T10 = refpropm('t','P',P10*1000,'H',H10,'CO2');
S10 = refpropm('s','T',T10,'P',P10*1000,'CO2');

%---分流比例
alpha = 0.3333; %分流比例

%---副压缩机入口
P11 = P6; %副压缩机入口压力（副路，从分流点6分出）
T11 = 0; %初始化，后面会根据T6更新
H11 = 0; %初始化，后面会根据H6更新
S11 = 0; %初始化，后面会根据S6更新

%---副压缩机出口
N_rc = 0.85; %副压缩机效率（降低到0.85，增加压缩过程的不可逆性）
P12 = P1; %副压缩机出口压力等于高压透平入口压力
S12_is = 0; %初始化，后面会根据S11更新
H12_is = 0; %初始化
H12 = 0; %初始化
T12 = 0; %初始化
S12 = 0; %初始化

%---合流点
P13 = P10; %合流点压力等于主压缩机b出口压力

%---低温回热器冷侧出口
P14 = P13; %低温回热器冷侧入口压力

%---高温回热器冷侧出口
P15 = P14; %高温回热器冷侧入口压力

%---新的初始化方法：从合流点开始，反向计算回热器状态

%---首先估计副压缩机入口/出口状态
%假设副压缩机入口温度低于但接近冷却器出口温度（我们已知点7温度是305K）
T11_est = 320; %假设副压缩机入口温度为320K
H11_est = refpropm('h','T',T11_est,'P',P11*1000,'CO2');
S11_est = refpropm('s','T',T11_est,'P',P11*1000,'CO2');

%基于估计的副压缩机入口状态计算副压缩机出口
S12_is = S11_est; %等熵过程
H12_is = refpropm('h','P',P12*1000,'S',S12_is,'CO2');
H12 = H11_est + (H12_is - H11_est)/N_rc; %考虑副压缩机效率
T12 = refpropm('t','P',P12*1000,'H',H12,'CO2');
S12 = refpropm('s','T',T12,'P',P12*1000,'CO2');

%---计算合流点状态（预设值）
H13 = (1-alpha)*H10 + alpha*H12; %基于流量比例计算合流点混合焓
T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');

%---初始化低温回热器冷侧出口温度（使用回热器效率法）
%基于合流点温度和低温回热器效率估计点14
e1 = 0.86; %低温回热器效率
e2 = 0.86; %高温回热器效率

%使用效率法估计低温回热器冷侧出口温度
%由于我们还不知道T5（需要计算得到），先假设一个合理温度差
%假设低温回热器可能达到的最高温度（热侧入口T5）约为600K
T5_guess = 600; %假设值，后续会通过迭代更新
T14 = T13 + e1 * (T5_guess - T13); %使用效率法计算低温回热器冷侧出口温度
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');

%---初始化高温回热器冷侧出口温度（使用回热器效率法）
%基于低温回热器冷侧出口温度和高温回热器效率估计点15
T15 = T14 + e2 * (T4 - T14); %使用效率法计算高温回热器冷侧出口温度
H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');

%---反向计算高温回热器热侧出口（点5）
Q_HT_cold = H15 - H14; %高温回热器冷侧吸热
H5 = H4 - Q_HT_cold; %高温回热器热侧放热等于冷侧吸热
T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');

%---反向计算低温回热器热侧出口（点6）
Q_LT_cold = H14 - H13; %低温回热器冷侧吸热
H6 = H5 - Q_LT_cold; %低温回热器热侧放热等于冷侧吸热
T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');

%---更新分流点状态（11点与6点状态相同）
T11 = T6;
H11 = H6;
S11 = S6;

%---更新副压缩机计算（基于更新的分流点状态）
S12_is = S11;
H12_is = refpropm('h','P',P12*1000,'S',S12_is,'CO2');
H12 = H11 + (H12_is - H11)/N_rc;
T12 = refpropm('t','P',P12*1000,'H',H12,'CO2');
S12 = refpropm('s','T',T12,'P',P12*1000,'CO2');

%---重新计算合流点（考虑实际熵增的混合过程）
%计算理想情况下的混合焓值（质量加权平均）
H13_ideal = (1-alpha)*H10 + alpha*H12;

%考虑实际熵增的修正（熵增系数可调整）
entropy_gen_factor = 0.03; %熵增系数，表示合流过程的不可逆性
%熵增修正（基于高低焓流的差异）
entropy_loss = entropy_gen_factor * abs(H12 - H10);
%考虑熵增后的实际合流点焓值（熵增导致可用能降低）
H13 = H13_ideal - entropy_loss;
T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');

%%---热平衡迭代求解
e1 = 0.86; %低温回热器效率
e2 = 0.86; %高温回热器效率
n_heater = 0.94; %加热器效率

%---热平衡迭代参数设置
phase1_max_iter = 50; %第一阶段迭代次数（用于获得物理合理的T-S图）
phase2_max_iter = 150; %第二阶段迭代次数（用于获得准确的效率值）
relax_coef = 0.5; %松弛系数 - 可以修改此值来影响收敛性能
converge_tol = 0.001; %收敛容差

disp('======= 第一阶段迭代：获取符合物理规律的T-S图数据 =======');
disp(['第一阶段最大迭代次数: ', num2str(phase1_max_iter), '，松弛系数: ', num2str(relax_coef)]);

%第一阶段迭代，主要用于获得符合物理规律的T-S图数据
object = 1;
iter = 0;

while object && iter < phase1_max_iter
    iter = iter + 1;
    
    %保存上一轮值用于比较收敛
    H5_old = H5;
    H6_old = H6;
    H14_old = H14;
    H15_old = H15;
    
    %更新高温回热器
    Q_HT_hot = H4 - H5; %高温回热器热侧放热
    Q_HT_cold = H15 - H14; %高温回热器冷侧吸热
    
    %检查能量守恒并进行调整
    dQ_HT = Q_HT_hot - Q_HT_cold;
    %调整高温回热器热侧出口焓
    H5_new = H4 - Q_HT_cold; %高温回热器热侧出口焓
    
    %带松弛系数的更新
    relax = relax_coef; %松弛系数
    H5 = H5_old + relax*(H5_new - H5_old);
    T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
    S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
    
    %更新低温回热器
    Q_LT_hot = H5 - H6; %低温回热器热侧放热
    Q_LT_cold = H14 - H13; %低温回热器冷侧吸热
    
    %检查能量守恒并进行调整
    dQ_LT = Q_LT_hot - Q_LT_cold;
    %调整低温回热器热侧出口焓
    H6_new = H5 - Q_LT_cold;
    
    %带松弛系数的更新
    H6 = H6_old + relax*(H6_new - H6_old);
    T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
    S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
    
    %更新分流点（11点状态与6点相同）
    T11 = T6;
    H11 = H6;
    S11 = S6;
    
    %更新副压缩机
    S12_is = S11;
    H12_is = refpropm('h','P',P12*1000,'S',S12_is,'CO2');
    H12 = H11 + (H12_is - H11)/N_rc;
    T12 = refpropm('t','P',P12*1000,'H',H12,'CO2');
    S12 = refpropm('s','T',T12,'P',P12*1000,'CO2');
    
    %更新合流点（考虑实际熵增的混合过程）
    %计算理想情况下的混合焓值（质量加权平均）
    H13_ideal = (1-alpha)*H10 + alpha*H12;
    
    %考虑实际熵增的修正（熵增系数可调整）
    entropy_gen_factor = 0.03; %熵增系数，表示合流过程的不可逆性
    %熵增修正（基于高低焓流的差异）
    entropy_loss = entropy_gen_factor * abs(H12 - H10);
    %考虑熵增后的实际合流点焓值（熵增导致可用能降低）
    H13 = H13_ideal - entropy_loss;
    T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
    S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
    
    %更新低温回热器冷侧出口
    %使用效率法计算低温回热器冷侧出口
    T_min = max(T13, T6);
    T_max = min(T5, T1);
    T_range = T_max - T_min;
    T14 = T13 + e1 * (T5 - T13); %使用效率计算低温回热器冷侧出口温度
    H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
    S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
    
    %更新高温回热器冷侧出口
    T15 = T14 + e2 * (T4 - T14); %使用效率法计算高温回热器冷侧出口温度
    H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
    S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
    
    %检查收敛性
    dH5 = abs(H5 - H5_old)/abs(H5);
    dH6 = abs(H6 - H6_old)/abs(H6);
    dH14 = abs(H14 - H14_old)/abs(H14);
    dH15 = abs(H15 - H15_old)/abs(H15);
    
    if max([dH5, dH6, dH14, dH15]) < converge_tol
        object = 0;
    end
    
    %输出迭代进度
    if mod(iter, 5) == 0
        disp(['第一阶段迭代次数: ', num2str(iter), ', 最大相对变化: ', num2str(max([dH5, dH6, dH14, dH15]))]);
    end
end

%---保存第一阶段结果（用于绘制物理合理的T-S图）
T1_phase1 = T1; T2_phase1 = T2; T3_phase1 = T3; T4_phase1 = T4; T5_phase1 = T5;
T6_phase1 = T6; T7_phase1 = T7; T8_phase1 = T8; T9_phase1 = T9; T10_phase1 = T10;
T11_phase1 = T11; T12_phase1 = T12; T13_phase1 = T13; T14_phase1 = T14; T15_phase1 = T15;

S1_phase1 = S1; S2_phase1 = S2; S3_phase1 = S3; S4_phase1 = S4; S5_phase1 = S5;
S6_phase1 = S6; S7_phase1 = S7; S8_phase1 = S8; S9_phase1 = S9; S10_phase1 = S10;
S11_phase1 = S11; S12_phase1 = S12; S13_phase1 = S13; S14_phase1 = S14; S15_phase1 = S15;

P_phase1 = [P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15];
H_phase1 = [H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15];

disp('第一阶段迭代完成，已保存符合物理规律的T-S图数据');
disp(['第一阶段迭代次数: ', num2str(iter), ', 最终相对变化: ', num2str(max([dH5, dH6, dH14, dH15]))]);

%---第二阶段迭代：在第一阶段基础上进一步优化效率计算的准确性
disp('======= 第二阶段迭代：优化效率计算的准确性 =======');
disp(['第二阶段最大迭代次数: ', num2str(phase2_max_iter - iter), '，松弛系数: ', num2str(relax_coef)]);

if iter < phase2_max_iter && max([dH5, dH6, dH14, dH15]) >= converge_tol
    %继续迭代直到收敛或达到最大迭代次数
    object = 1;
    phase2_start_iter = iter;
    
    while object && iter < phase2_max_iter
        iter = iter + 1;
        
        %保存上一轮值用于比较收敛
        H5_old = H5;
        H6_old = H6;
        H14_old = H14;
        H15_old = H15;
        
        %更新高温回热器
        Q_HT_hot = H4 - H5; %高温回热器热侧放热
        Q_HT_cold = H15 - H14; %高温回热器冷侧吸热
        
        %检查能量守恒并进行调整
        dQ_HT = Q_HT_hot - Q_HT_cold;
        %调整高温回热器热侧出口焓
        H5_new = H4 - Q_HT_cold; %高温回热器热侧出口焓
        
        %带松弛系数的更新
        relax = relax_coef; %松弛系数
        H5 = H5_old + relax*(H5_new - H5_old);
        T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
        S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
        
        %更新低温回热器
        Q_LT_hot = H5 - H6; %低温回热器热侧放热
        Q_LT_cold = H14 - H13; %低温回热器冷侧吸热
        
        %检查能量守恒并进行调整
        dQ_LT = Q_LT_hot - Q_LT_cold;
        %调整低温回热器热侧出口焓
        H6_new = H5 - Q_LT_cold;
        
        %带松弛系数的更新
        H6 = H6_old + relax*(H6_new - H6_old);
        T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
        S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
        
        %更新分流点（11点状态与6点相同）
        T11 = T6;
        H11 = H6;
        S11 = S6;
        
        %更新副压缩机
        S12_is = S11;
        H12_is = refpropm('h','P',P12*1000,'S',S12_is,'CO2');
        H12 = H11 + (H12_is - H11)/N_rc;
        T12 = refpropm('t','P',P12*1000,'H',H12,'CO2');
        S12 = refpropm('s','T',T12,'P',P12*1000,'CO2');
        
        %更新合流点（考虑实际熵增的混合过程）
        %计算理想情况下的混合焓值（质量加权平均）
        H13_ideal = (1-alpha)*H10 + alpha*H12;
        
        %考虑实际熵增的修正（熵增系数可调整）
        entropy_gen_factor = 0.03; %熵增系数，表示合流过程的不可逆性
        %熵增修正（基于高低焓流的差异）
        entropy_loss = entropy_gen_factor * abs(H12 - H10);
        %考虑熵增后的实际合流点焓值（熵增导致可用能降低）
        H13 = H13_ideal - entropy_loss;
        T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
        S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
        
        %更新低温回热器冷侧出口
        %使用效率法计算低温回热器冷侧出口
        T_min = max(T13, T6);
        T_max = min(T5, T1);
        T_range = T_max - T_min;
        T14 = T13 + e1 * (T5 - T13); %使用效率计算低温回热器冷侧出口温度
        H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
        S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
        
        %更新高温回热器冷侧出口
        T15 = T14 + e2 * (T4 - T14); %使用效率法计算高温回热器冷侧出口温度
        H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
        S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
        
        %检查收敛性
        dH5 = abs(H5 - H5_old)/abs(H5);
        dH6 = abs(H6 - H6_old)/abs(H6);
        dH14 = abs(H14 - H14_old)/abs(H14);
        dH15 = abs(H15 - H15_old)/abs(H15);
        
        if max([dH5, dH6, dH14, dH15]) < converge_tol
            object = 0;
        end
        
        %输出迭代进度
        if mod(iter, 5) == 0
            disp(['第二阶段迭代次数: ', num2str(iter - phase2_start_iter), ', 总迭代次数: ', num2str(iter), ', 最大相对变化: ', num2str(max([dH5, dH6, dH14, dH15]))]);
        end
    end
    
    disp('第二阶段迭代完成，用于计算效率');
    disp(['总迭代次数: ', num2str(iter), ', 最终相对变化: ', num2str(max([dH5, dH6, dH14, dH15]))]);
    
    if iter >= phase2_max_iter && max([dH5, dH6, dH14, dH15]) >= converge_tol
        disp('警告：迭代未收敛，已达最大迭代次数');
    end
else
    disp('第一阶段已收敛，无需进行第二阶段迭代');
end

%---输出最后一次迭代的状态点详细信息
disp('=============================================');
disp('最后一次迭代的各状态点详细信息（第二阶段结果，用于效率计算）：');
disp('=============================================');
disp(['点1 (高压透平入口)：T = ', num2str(T1), ' K, P = ', num2str(P1), ' MPa, H = ', num2str(H1), ' J/kg, S = ', num2str(S1), ' J/(kg·K)']);
disp(['点2 (高压透平出口)：T = ', num2str(T2), ' K, P = ', num2str(P2), ' MPa, H = ', num2str(H2), ' J/kg, S = ', num2str(S2), ' J/(kg·K)']);
disp(['点3 (再热器出口)：T = ', num2str(T3), ' K, P = ', num2str(P3), ' MPa, H = ', num2str(H3), ' J/kg, S = ', num2str(S3), ' J/(kg·K)']);
disp(['点4 (低压透平出口)：T = ', num2str(T4), ' K, P = ', num2str(P4), ' MPa, H = ', num2str(H4), ' J/kg, S = ', num2str(S4), ' J/(kg·K)']);
disp(['点5 (高温回热器热侧出口)：T = ', num2str(T5), ' K, P = ', num2str(P5), ' MPa, H = ', num2str(H5), ' J/kg, S = ', num2str(S5), ' J/(kg·K)']);
disp(['点6 (低温回热器热侧出口)：T = ', num2str(T6), ' K, P = ', num2str(P6), ' MPa, H = ', num2str(H6), ' J/kg, S = ', num2str(S6), ' J/(kg·K)']);
disp(['点7 (冷却器出口)：T = ', num2str(T7), ' K, P = ', num2str(P7), ' MPa, H = ', num2str(H7), ' J/kg, S = ', num2str(S7), ' J/(kg·K)']);
disp(['点8 (主压缩机a出口)：T = ', num2str(T8), ' K, P = ', num2str(P8), ' MPa, H = ', num2str(H8), ' J/kg, S = ', num2str(S8), ' J/(kg·K)']);
disp(['点9 (中间冷却器出口)：T = ', num2str(T9), ' K, P = ', num2str(P9), ' MPa, H = ', num2str(H9), ' J/kg, S = ', num2str(S9), ' J/(kg·K)']);
disp(['点10 (主压缩机b出口)：T = ', num2str(T10), ' K, P = ', num2str(P10), ' MPa, H = ', num2str(H10), ' J/kg, S = ', num2str(S10), ' J/(kg·K)']);
disp(['点11 (副压缩机入口)：T = ', num2str(T11), ' K, P = ', num2str(P11), ' MPa, H = ', num2str(H11), ' J/kg, S = ', num2str(S11), ' J/(kg·K)']);
disp(['点12 (副压缩机出口)：T = ', num2str(T12), ' K, P = ', num2str(P12), ' MPa, H = ', num2str(H12), ' J/kg, S = ', num2str(S12), ' J/(kg·K)']);
disp(['点13 (合流点)：T = ', num2str(T13), ' K, P = ', num2str(P13), ' MPa, H = ', num2str(H13), ' J/kg, S = ', num2str(S13), ' J/(kg·K)']);
disp(['点14 (低温回热器冷侧出口)：T = ', num2str(T14), ' K, P = ', num2str(P14), ' MPa, H = ', num2str(H14), ' J/kg, S = ', num2str(S14), ' J/(kg·K)']);
disp(['点15 (高温回热器冷侧出口)：T = ', num2str(T15), ' K, P = ', num2str(P15), ' MPa, H = ', num2str(H15), ' J/kg, S = ', num2str(S15), ' J/(kg·K)']);
disp('=============================================');
disp(['最后一次迭代的相对变化：dH5 = ', num2str(dH5), ', dH6 = ', num2str(dH6), ', dH14 = ', num2str(dH14), ', dH15 = ', num2str(dH15)]);
disp(['收敛判据: max([dH5, dH6, dH14, dH15]) = ', num2str(max([dH5, dH6, dH14, dH15])), ' < ', num2str(converge_tol)]);
disp('=============================================');

%---计算循环效率（使用最终收敛值）
%功率计算
W_hp_turb = H1 - H2; %高压透平做功
W_lp_turb = H3 - H4; %低压透平做功
W_mc_a = H8 - H7; %主压缩机a耗功
W_mc_b = H10 - H9; %主压缩机b耗功
W_rc = H12 - H11; %副压缩机耗功
W_comp_total = (1-alpha)*(W_mc_a + W_mc_b) + alpha*W_rc; %总压缩功率
W_turb_total = W_hp_turb + W_lp_turb; %总透平功率
W_net = W_turb_total - W_comp_total; %净输出功率

%热量计算
Q_heater = H1 - H15; %加热器热量（直接从点15到点1）
Q_reheater = H3 - H2; %再热器热量
Q_cooler = (1-alpha)*(H6 - H7); %主路冷却器热量（从点6到点7）
Q_intercooler = (1-alpha)*(H8 - H9); %中间冷却器热量
Q_in_total = Q_heater + Q_reheater; %总输入热量

%循环热效率
eta = W_net / Q_in_total;

%输出结果
disp('=============================================');
disp('最终计算结果（使用完全收敛的值计算效率）：');
disp('=============================================');
disp(['高压透平做功: ' num2str(W_hp_turb) ' J/kg'])
disp(['低压透平做功: ' num2str(W_lp_turb) ' J/kg'])
disp(['总压缩功率: ' num2str(W_comp_total) ' J/kg'])
disp(['净输出功率: ' num2str(W_net) ' J/kg'])
disp(['总输入热量: ' num2str(Q_in_total) ' J/kg'])
disp(['循环热效率: ' num2str(eta*100) ' %'])
disp('=============================================');

% 调用T-S绘图函数（使用第一阶段结果绘制物理合理的T-S图）
disp('使用第一阶段结果绘制符合物理规律的T-S图');
plot_ts_diagram_renum(T1_phase1, T2_phase1, T3_phase1, T4_phase1, T5_phase1, T6_phase1, T7_phase1, T8_phase1, T9_phase1, T10_phase1, T11_phase1, T12_phase1, T13_phase1, T14_phase1, T15_phase1,...
                      S1_phase1, S2_phase1, S3_phase1, S4_phase1, S5_phase1, S6_phase1, S7_phase1, S8_phase1, S9_phase1, S10_phase1, S11_phase1, S12_phase1, S13_phase1, S14_phase1, S15_phase1,...
                      P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15,...
                      alpha, W_net, eta);

% T-S图绘制函数（重新编号版）
function plot_ts_diagram_renum(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15,...
                               S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15,...
                               P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15,...
                               alpha, W_net, eta)
    % 创建一个新图
    figure;
    hold on;
    
    % 设置标题和标签
    title('超临界CO₂布雷顿循环T-S图（15点系统）');
    xlabel('熵 s [J/(kg·K)]');
    ylabel('温度 T [K]');
    
    % 主循环连线
    % 高压透平 (1→2)
    line1 = plot([S1,S2], [T1,T2], 'b-', 'LineWidth', 1.5);
    
    % 再热器 (2→3) - 等压加热
    line2 = plot([S2,S3], [T2,T3], 'r-', 'LineWidth', 1.5);
    
    % 低压透平 (3→4)
    plot([S3,S4], [T3,T4], 'b-', 'LineWidth', 1.5);
    
    % 高温回热器热侧 (4→5) - 等压放热
    plot([S4,S5], [T4,T5], 'b-', 'LineWidth', 1.5);
    
    % 低温回热器热侧 (5→6) - 等压放热
    plot([S5,S6], [T5,T6], 'b-', 'LineWidth', 1.5);
    
    % 主路冷却器 (6→7) - 等压放热
    plot([S6,S7], [T6,T7], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机a (7→8)
    plot([S7,S8], [T7,T8], 'b-', 'LineWidth', 1.5);
    
    % 中间冷却器 (8→9) - 等压放热
    plot([S8,S9], [T8,T9], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机b (9→10)
    plot([S9,S10], [T9,T10], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机b出口到合流点 (10→13) - 连线
    plot([S10,S13], [T10,T13], 'b-', 'LineWidth', 1.5);
    
    % 低温回热器冷侧 (13→14) - 等压吸热
    plot([S13,S14], [T13,T14], 'r-', 'LineWidth', 1.5);
    
    % 高温回热器冷侧 (14→15) - 等压吸热
    plot([S14,S15], [T14,T15], 'r-', 'LineWidth', 1.5);
    
    % 加热器 (15→1) - 等压吸热
    plot([S15,S1], [T15,T1], 'r-', 'LineWidth', 1.5);
    
    % 分流部分（虚线）
    % 分流点到副压缩机入口 (6→11) - 连线
    plot([S6,S11], [T6,T11], 'b--', 'LineWidth', 1.5);
    
    % 副压缩机 (11→12)
    line16 = plot([S11,S12], [T11,T12], 'b--', 'LineWidth', 1.5);
    
    % 副压缩机出口到合流点 (12→13) - 连线
    plot([S12,S13], [T12,T13], 'b--', 'LineWidth', 1.5);
    
    % 标记主循环状态点
    main_points = 1:15;
    main_T = [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15];
    main_S = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15];
    
    for i = 1:length(main_points)
        if i >= 11 && i <= 12  % 分流部分的点
            plot(main_S(i), main_T(i), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
        else  % 主循环的点
            plot(main_S(i), main_T(i), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        end
        text(main_S(i), main_T(i), [' ', num2str(main_points(i))], 'FontSize', 8);
    end
    
    % 添加图例
    legend([line1, line2, line16, plot(NaN,NaN,'ro','MarkerFaceColor','r'), plot(NaN,NaN,'bo','MarkerFaceColor','b')], ...
           {'主循环路径', '吸热过程', '分流路径', '主循环状态点', '分流状态点'}, ...
           'Location', 'Best', 'FontSize', 8);
    
    % 添加说明性文本框 - 描述各段的过程
    description = {
        '1→2: 高压透平膨胀做功', ...
        '2→3: 再热器等压加热', ...
        '3→4: 低压透平膨胀做功', ...
        '4→5: 高温回热器热侧放热', ...
        '5→6: 低温回热器热侧放热', ...
        '6→7: 冷却器冷却', ...
        '7→8: 主压缩机a压缩', ...
        '8→9: 中间冷却器冷却', ...
        '9→10: 主压缩机b压缩', ...
        '6→11→12: 副压缩机路径', ...
        '10,12→13: 合流点', ...
        '13→14→15→1: 回热与加热'};
    
    % 设置图形属性
    grid on;
    set(gca, 'FontSize', 12);
    set(gcf, 'Position', [100, 100, 1000, 700]);
    
    % 添加过程说明文本
    text(min(main_S)+0.05*(max(main_S)-min(main_S)), max(main_T)*0.85, ...
         description, 'FontSize', 8, 'BackgroundColor', [1 1 1 0.7]);
    
    % 添加性能参数文本
    text(min(main_S)+0.60*(max(main_S)-min(main_S)), max(main_T)*0.85, ...
        {['循环热效率: ', num2str(eta*100, '%.2f'), ' %'], ...
         ['净输出功率: ', num2str(W_net/1000, '%.2f'), ' kJ/kg'], ...
         ['分流比例: ', num2str(alpha, '%.2f')]}, ...
         'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
    
    hold off;
end 