%%------初始参数
%%---透平
n2 = 0.93;%透平效率
P14 = str2double('21');   %进口压力21，单位Mpa（随模型变化）
T14 = str2double('873');  %进口温度，873K;600℃（随模型变化）
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
P1 = str2double('7.5729');%出口压力，单位Mpa（随模型变化）
S1 = S14;
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');%点1的理想焓
T1 = refpropm('t','H',H1,'P',P1*1000,'CO2');%点1的温度
H1 = H14 + (H1 - H14)*n2;                   %点1实际的焓
S1 = refpropm('s','P',P1*1000,'H',H1,'CO2');%点1实际的熵

%%---主压缩机（1-x分流，效率N1）
N1 = 0.89;
T5 = str2double('305'); %主压缩机进口温度32℃（随模型变化）
P5 = P1;                %主压缩机进口压力，P5等于透平出口压力P1
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
H5 = refpropm('h','T',T5,'P',P5*1000,'CO2');
S6 = S5;                %假设压缩完熵不变
P6 = P14;               %主压缩机出口与透平进口相连，压力相同
H6 = refpropm('h','P',P6*1000,'S',S6,'CO2');
H6 = H5 + (H6 - H5)/N1;%点6实际焓值
T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');%点6实际温度
S6 = refpropm('s','T',T6,'H',H6,'CO2');%点6实际的熵

%%---副压缩机（x分流,效率N2）
N2 = 0.89;
P9 = P1;                 %副压缩机进口与透平出口相连，压力相等

%%---冷却器（1-x分流）
P4 = P1;%冷却器进口与透平出口相连，压力相等

%%低温回热器
T7 = T6;
P7 = P6;
H7 = H6;
P8 = P7;

%%高温回热器
P10 = P14;
P11 = P10;
P3 = P1;
T2 = T1;
H2 = refpropm('h','T',T2,'P',P1*1000,'CO2');
P12 = P11;

%%---求解
e1 = 0.86;%低温回热器
e2 = 0.86;%高温回热器
n1 = 0.94;%锅炉效率
H111 = 300000;%初始假设点11的焓值

THTR = 15;%高温回热器冷端端差
TLTR = 15;%低温回热器冷端端差

x = 0.3333;
    object = 1;
    while object
        T11 = refpropm('t','H',H111,'P',P11*1000,'CO2');
        T3 = T11 + THTR;%高温回热器热端出口（利用端差计算）
        H3 = refpropm('H','T',T3,'P',P3*1000,'CO2');
        H12 = H2 - H3 + H111;
        T12 = refpropm('t','H',H12,'P',P12*1000,'CO2');
        T3 = refpropm('t','H',H3,'P',P3*1000,'CO2');
        T8 = T3 - TLTR;
        H8 = refpropm('h','T',T8,'P',P8*1000,'CO2');
        H9 = H3 - (1-x)*(H8 - H7);
        T9 = refpropm('t','H',H9,'P',P9*1000,'CO2');
        S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');
        S10 = S9;%        
        P10 = P14;
        H10 = refpropm('h','P',P10*1000,'S',S10,'CO2');
        H10 = H9 + (H10 - H9)/N2;%点10实际焓值
        T10 = refpropm('t','P',P10*1000,'S',S10,'CO2');%点10实际温度
        S10 = refpropm('s','T',T10,'H',H10,'CO2');%点10实际熵
        %冷却器
        H4 = H9;
        %分流汇合
        H11 = x*H10 + (1-x)*H8;
        if (0.999 <= H11/H111) && (H11/H111<= 1.001)
            object = 0;
        else
            H111 = H111 + 50;
        end
    end
n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1); %循环效率
%-----环境温度25℃
disp(['循环效率',num2str(n)])

%% 绘制T-S图
%> @brief 绘制超临界二氧化碳循环的T-S图
%> @author Claude AI
%> @date 自动生成

% 收集所有状态点的温度和熵值
T_points = [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T14];
S_points = [S1, refpropm('s','T',T2,'P',P3*1000,'CO2'), refpropm('s','T',T3,'P',P3*1000,'CO2'), ...
           refpropm('s','T',T4,'P',P4*1000,'CO2'), S5, S6, ...
           refpropm('s','T',T7,'P',P7*1000,'CO2'), refpropm('s','T',T8,'P',P8*1000,'CO2'), ...
           S9, S10, refpropm('s','T',T11,'P',P11*1000,'CO2'), ...
           refpropm('s','T',T12,'P',P12*1000,'CO2'), S14];
point_labels = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '14'};

% 创建图形窗口
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% 绘制状态点
scatter(S_points, T_points, 70, 'ro', 'filled');
hold on;

% 添加状态点标签
for i = 1:length(T_points)
    text(S_points(i)+0.01, T_points(i)+5, point_labels{i}, 'FontSize', 12);
end

% 连接状态点，显示热力循环过程
% 主要过程
plot([S14, S1], [T14, T1], 'b-', 'LineWidth', 2); % 透平膨胀过程 14->1
plot([S5, S6], [T5, T6], 'b-', 'LineWidth', 2);   % 主压缩机压缩过程 5->6
plot([S9, S10], [T9, T10], 'b-', 'LineWidth', 2); % 副压缩机压缩过程 9->10

% 热交换过程
S_htr_hot = linspace(S_points(2), S_points(3), 20);
T_htr_hot = zeros(size(S_htr_hot));
for i = 1:length(S_htr_hot)
    T_htr_hot(i) = refpropm('t', 'p', P3*1000, 's', S_htr_hot(i), 'CO2');
end
plot(S_htr_hot, T_htr_hot, 'r-', 'LineWidth', 2); % 高温回热器热侧 2->3

S_htr_cold = linspace(S_points(11), S_points(12), 20);
T_htr_cold = zeros(size(S_htr_cold));
for i = 1:length(S_htr_cold)
    T_htr_cold(i) = refpropm('t', 'p', P11*1000, 's', S_htr_cold(i), 'CO2');
end
plot(S_htr_cold, T_htr_cold, 'r-', 'LineWidth', 2); % 高温回热器冷侧 11->12

% 绘制恒压线作为参考
% 高压侧
S_highP = linspace(min(S_points)*0.95, max(S_points)*1.05, 100);
T_highP = zeros(size(S_highP));
for i = 1:length(S_highP)
    try
        T_highP(i) = refpropm('t', 'p', P14*1000, 's', S_highP(i), 'CO2');
    catch
        T_highP(i) = NaN;
    end
end
plot(S_highP, T_highP, 'k--', 'LineWidth', 1);

% 低压侧
S_lowP = linspace(min(S_points)*0.95, max(S_points)*1.05, 100);
T_lowP = zeros(size(S_lowP));
for i = 1:length(S_lowP)
    try
        T_lowP(i) = refpropm('t', 'p', P1*1000, 's', S_lowP(i), 'CO2');
    catch
        T_lowP(i) = NaN;
    end
end
plot(S_lowP, T_lowP, 'k--', 'LineWidth', 1);

% 添加图形标题和标签
title('超临界二氧化碳布雷顿循环T-S图', 'FontSize', 16);
xlabel('熵 S (J/kg·K)', 'FontSize', 14);
ylabel('温度 T (K)', 'FontSize', 14);
grid on;

% 添加图例
legend('状态点', '透平膨胀', '主压缩机', '副压缩机', '高温回热器(热)', '高温回热器(冷)', ...
       [num2str(P14), ' MPa'], [num2str(P1), ' MPa'], 'Location', 'best');

% 显示循环效率
efficiency_text = sprintf('循环效率: %.2f%%', n*100);
text(min(S_points), max(T_points)*0.9, efficiency_text, 'FontSize', 14, 'BackgroundColor', [0.9, 0.9, 0.9]);

hold off;

