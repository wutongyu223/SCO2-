%%------初始参数
%---高压透平
n_hp_turb = 0.93; %高压透平效率
P1 = str2double('21'); %高压透平进口压力21，单位MPa
T1 = str2double('873'); %高压透平进口温度，873K;600℃
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');
H1 = refpropm("H","T",T1,"P",P1*1000,"CO2")

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

%---低压透平
n_lp_turb = 0.93; %低压透平效率
P4 = str2double('7.5729'); %低压透平出口压力，单位MPa
S4_is = S3; %等熵过程
H4_is = refpropm('h','P',P4*1000,'S',S4_is,'CO2'); %等熵焓变
H4 = H3 - (H3 - H4_is)*n_lp_turb; %实际焓变
T4 = refpropm('t','P',P4*1000,'H',H4,'CO2'); %点4的温度
S4 = refpropm('s','T',T4,'P',P4*1000,'CO2'); %点4实际的熵

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
    
    %更新低温回热器
    Q_LT_hot = H5 - H6; %低温回热器热侧放热
    Q_LT_cold = H15 - H14; %低温回热器冷侧吸热
    
    %检查能量守恒并进行调整
    dQ_LT = Q_LT_hot - Q_LT_cold;
    %调整低温回热器热侧出口焓
    H6_new = H5 - Q_LT_cold;
    
    %带松弛系数的更新
    H6 = H6_old + relax*(H6_new - H6_old);
    T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
    S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
    
    %更新分流点
    T7 = T6;
    T12 = T6;
    H7 = H6;
    H12 = H6;
    S7 = S6;
    S12 = S6;
    
    %更新副压缩机
    S13_is = S12;
    H13_is = refpropm('h','P',P13*1000,'S',S13_is,'CO2');
    H13 = H12 + (H13_is - H12)/N_rc;
    T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
    S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
    
    %更新合流点
    H14 = (1-alpha)*H11 + alpha*H13;
    T14 = refpropm('t','P',P14*1000,'H',H14,'CO2');
    S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
    
    %更新低温回热器冷侧出口
    %使用效率法计算低温回热器冷侧出口
    T_min = max(T14, T6);
    T_max = min(T5, T17);
    T_range = T_max - T_min;
    T15 = T14 + e1 * (T5 - T14); %使用效率计算低温回热器冷侧出口温度
    H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
    S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
    
    %更新高温回热器冷侧出口
    T16 = T15 + e2 * (T4 - T15); %使用效率计算高温回热器冷侧出口温度
    H16 = refpropm('h','T',T16,'P',P16*1000,'CO2');
    S16 = refpropm('s','T',T16,'P',P16*1000,'CO2');
    
    %检查收敛性
    dH5 = abs(H5 - H5_old)/abs(H5);
    dH6 = abs(H6 - H6_old)/abs(H6);
    dH15 = abs(H15 - H15_old)/abs(H15);
    dH16 = abs(H16 - H16_old)/abs(H16);
    
    if max([dH5, dH6, dH15, dH16]) < 0.001
        object = 0;
    end
    
    %输出迭代进度
    if mod(iter, 5) == 0
        disp(['迭代次数: ', num2str(iter), ', 最大相对变化: ', num2str(max([dH5, dH6, dH15, dH16]))]);
    end
end

if iter >= max_iter
    disp('警告：迭代未收敛，已达最大迭代次数');
end

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

% T-S图绘制函数放在文件末尾
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
    
    % 低压透平 (3→4)
    line3 = plot([S3,S4], [T3,T4], 'b-', 'LineWidth', 1.5);
    
    % 高温回热器热侧 (4→5) - 等压放热
    line4 = plot([S4,S5], [T4,T5], 'b-', 'LineWidth', 1.5);
    
    % 低温回热器热侧 (5→6) - 等压放热
    line5 = plot([S5,S6], [T5,T6], 'b-', 'LineWidth', 1.5);
    
    % 分流点到主路冷却器 (6→7) - 连线
    line6 = plot([S6,S7], [T6,T7], 'b-', 'LineWidth', 1.5);
    
    % 主路冷却器 (7→8) - 等压放热
    line7 = plot([S7,S8], [T7,T8], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机a (8→9)
    line8 = plot([S8,S9], [T8,T9], 'b-', 'LineWidth', 1.5);
    
    % 中间冷却器 (9→10) - 等压放热
    line9 = plot([S9,S10], [T9,T10], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机b (10→11)
    line10 = plot([S10,S11], [T10,T11], 'b-', 'LineWidth', 1.5);
    
    % 主压缩机b出口到合流点 (11→14) - 连线
    line11 = plot([S11,S14], [T11,T14], 'b-', 'LineWidth', 1.5);
    
    % 低温回热器冷侧 (14→15) - 等压吸热
    line12 = plot([S14,S15], [T14,T15], 'r-', 'LineWidth', 1.5);
    
    % 高温回热器冷侧 (15→16) - 等压吸热
    line13 = plot([S15,S16], [T15,T16], 'r-', 'LineWidth', 1.5);
    
    % 加热器 (16→17/1) - 等压吸热
    line14 = plot([S16,S17], [T16,T17], 'r-', 'LineWidth', 1.5);
    
    % 分流部分（虚线）
    % 分流点到副压缩机入口 (6→12) - 连线
    line15 = plot([S6,S12], [T6,T12], 'b--', 'LineWidth', 1.5);
    
    % 副压缩机 (12→13)
    line16 = plot([S12,S13], [T12,T13], 'b--', 'LineWidth', 1.5);
    
    % 副压缩机出口到合流点 (13→14) - 连线
    line17 = plot([S13,S14], [T13,T14], 'b--', 'LineWidth', 1.5);
    
    % 标记主循环状态点
    main_points = [1,2,3,4,5,6,7,8,9,10,11,14,15,16,17];
    main_T = [T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T14,T15,T16,T17];
    main_S = [S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S14,S15,S16,S17];
    
    for i = 1:length(main_points)
        plot(main_S(i), main_T(i), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        text(main_S(i), main_T(i), [' ', num2str(main_points(i))], 'FontSize', 8);
    end
    
    % 标记分流部分的状态点（12和13）
    plot(S12, T12, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    text(S12, T12, ' 12', 'FontSize', 8);
    plot(S13, T13, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    text(S13, T13, ' 13', 'FontSize', 8);
    
    % 添加图例
    legend([line1, line2, line16, plot(NaN,NaN,'ro','MarkerFaceColor','r'), plot(NaN,NaN,'bo','MarkerFaceColor','b')], ...
           {'主循环路径', '吸热过程', '分流路径', '主循环状态点', '分流状态点'}, ...
           'Location', 'Best', 'FontSize', 8);
    
    % 设置图形属性
    grid on;
    set(gca, 'FontSize', 12);
    set(gcf, 'Position', [100, 100, 800, 600]);
    
    % 添加说明文本
    text(min(main_S)+0.05*(max(main_S)-min(main_S)), max(main_T)*0.9, ...
        {['循环热效率: ', num2str(eta*100, '%.2f'), ' %'], ...
         ['净输出功率: ', num2str(W_net/1000, '%.2f'), ' kJ/kg']}, ...
         'FontSize', 10, 'BackgroundColor', [1 1 1 0.7]);
    
    hold off;
end

% 辅助函数：绘制等压过程线 (不再使用，采用简单直线)
function draw_isobaric_process(T1, S1, T2, S2, P, linespec)
    % 在两点之间绘制等压线，使用更多点进行平滑插值
    n_points = 20;
    
    % 如果两点温度非常接近，直接连线
    if abs(T1 - T2) < 1
        plot([S1, S2], [T1, T2], linespec);
        return;
    end
    
    % 否则创建中间点
    T_points = linspace(T1, T2, n_points);
    S_points = zeros(1, n_points);
    
    for i = 1:n_points
        % 获取给定温度和压力下的熵
        S_points(i) = refpropm('s', 'T', T_points(i), 'P', P*1000, 'CO2');
    end
    
    % 绘制等压线
    plot(S_points, T_points, linespec);
end

