%%------初始参数
%%---高压透平
n_hp_t = 0.93; %高压透平效率
P1 = str2double('21');   %高压透平入口压力，单位MPa
T1 = str2double('873');  %高压透平入口温度，873K;600℃
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');

%%---再热段
P2 = str2double('15');  %高压透平出口/再热器入口压力，单位MPa
S2 = S1;               %假设等熵膨胀
H2_is = refpropm('h','P',P2*1000,'S',S2,'CO2'); %理想等熵焓
H2 = H1 + (H2_is - H1)*n_hp_t;   %实际焓值
T2 = refpropm('t','P',P2*1000,'H',H2,'CO2');
S2 = refpropm('s','T',T2,'P',P2*1000,'CO2');  %实际熵

%%---再热器
P3 = P2;                %再热器出口/低压透平入口压力
T3 = str2double('873'); %再热温度，873K
S3 = refpropm('s','T',T3,'P',P3*1000,'CO2');
H3 = refpropm('h','T',T3,'P',P3*1000,'CO2');

%%---低压透平
n_lp_t = 0.93; %低压透平效率
P4 = str2double('7.6');  %低压透平出口/高温回热器热侧入口压力，单位MPa
S4 = S3;                %假设等熵膨胀
H4_is = refpropm('h','P',P4*1000,'S',S4,'CO2'); %理想等熵焓
H4 = H3 + (H4_is - H3)*n_lp_t;   %实际焓值
T4 = refpropm('t','P',P4*1000,'H',H4,'CO2');
S4 = refpropm('s','T',T4,'P',P4*1000,'CO2');  %实际熵

%%---高温回热器热侧
P5 = P4; %高温回热器热侧出口/低温回热器热侧入口压力
% 初始猜测T5值
T5 = 450; % 初始猜测值，将在迭代中更新
H5 = refpropm('h','T',T5,'P',P5*1000,'CO2');
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');

%%---低温回热器热侧
P6 = P5; %低温回热器热侧出口/分流点压力
% 初始猜测T6值
T6 = 350; % 初始猜测值，将在迭代中更新
H6 = refpropm('h','T',T6,'P',P6*1000,'CO2');
S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');

%%---主路冷却器
P7 = P6; %冷却器入口压力（主路）
T7 = T6; % 初始值等于T6
H7 = H6;
S7 = S6;

%%---主压缩机a
N_mc_a = 0.89; %主压缩机a效率
P8 = P7;      %冷却器出口/主压缩机a入口压力
T8 = str2double('305'); %冷却器出口温度，305K
S8 = refpropm('s','T',T8,'P',P8*1000,'CO2');
H8 = refpropm('h','T',T8,'P',P8*1000,'CO2');

%%---中间冷却器
P9 = str2double('15'); %主压缩机a出口/中间冷却入口压力
S9 = S8;               %假设等熵压缩
H9_is = refpropm('h','P',P9*1000,'S',S9,'CO2'); %理想等熵焓
H9 = H8 + (H9_is - H8)/N_mc_a;   %实际焓值
T9 = refpropm('t','P',P9*1000,'H',H9,'CO2');
S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');  %实际熵

%%---主压缩机b
N_mc_b = 0.89; %主压缩机b效率
P10 = P9;      %中间冷却出口/主压缩机b入口压力
T10 = str2double('320'); %中间冷却出口温度（中间冷却后的温度）
S10 = refpropm('s','T',T10,'P',P10*1000,'CO2');
H10 = refpropm('h','T',T10,'P',P10*1000,'CO2');
P11 = P1;      %主压缩机b出口/合流点压力
S11 = S10;     %假设等熵压缩
H11_is = refpropm('h','P',P11*1000,'S',S11,'CO2'); %理想等熵焓
H11 = H10 + (H11_is - H10)/N_mc_b;   %实际焓值
T11 = refpropm('t','P',P11*1000,'H',H11,'CO2');
S11 = refpropm('s','T',T11,'P',P11*1000,'CO2');  %实际熵

%%---副压缩机（分流路径）
N_rc = 0.89;   %副压缩机效率
P12 = P6;      %分流点/副压缩机入口压力
T12 = T6;      % 初始值等于T6
H12 = H6;
S12 = S6;
P13 = P1;      %副压缩机出口/合流点压力
% 初始计算副压缩机出口状态
S13_is = S12;  % 理想等熵过程
H13_is = refpropm('h','P',P13*1000,'S',S13_is,'CO2');
H13 = H12 + (H13_is - H12)/N_rc; % 实际焓值
T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');

%%---合流段
P14 = P11;     %合流点出口/低温回热器冷侧入口压力
% 初始计算合流点状态
alpha = 0.35;  %分流比例（分流到副压缩机的比例）
H14_init = alpha*H13 + (1-alpha)*H11;
T14 = refpropm('t','H',H14_init,'P',P14*1000,'CO2');
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');

%%---回热器参数
e_ltr = 0.86;  %低温回热器效率
e_htr = 0.86;  %高温回热器效率
n_heater = 0.94; %加热器效率

%%---迭代初始值
THTR = 15;     %高温回热器冷端端差
TLTR = 15;     %低温回热器冷端端差
H14_guess = H14_init; %初始猜测合流点焓值

%%---求解
object = 1;
iteration = 0;
max_iterations = 100;

while object && (iteration < max_iterations)
    iteration = iteration + 1;
    
    % 合流点计算
    H14 = H14_guess;
    T14 = refpropm('t','H',H14,'P',P14*1000,'CO2');
    S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
    
    % 低温回热器冷侧出口/高温回热器冷侧入口
    P15 = P14;
    % 使用低温回热器冷端端差
    T15 = T6 - TLTR; 
    H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
    S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
    
    % 高温回热器冷侧出口/加热器入口
    P16 = P15;
    % 使用高温回热器冷端端差
    T16 = T5 - THTR;
    H16 = refpropm('h','T',T16,'P',P16*1000,'CO2');
    S16 = refpropm('s','T',T16,'P',P16*1000,'CO2');
    
    % 加热器出口（等于高压透平入口）
    P17 = P1;
    T17 = T1;
    H17 = H1;
    S17 = S1;
    
    % 低温回热器传热量计算
    Q_ltr = H15 - H14; % 低温回热器传热量（冷侧）
    
    % 高温回热器传热量计算
    Q_htr = H16 - H15; % 高温回热器传热量（冷侧）
    
    % 高温回热器热侧计算
    H5_new = H4 - Q_htr;   % 高温回热器热侧出口焓
    T5_new = refpropm('t','H',H5_new,'P',P5*1000,'CO2');
    S5_new = refpropm('s','T',T5_new,'P',P5*1000,'CO2');
    
    % 低温回热器热侧计算
    H6_new = H5_new - Q_ltr/(1-alpha); % 低温回热器热侧出口焓
    T6_new = refpropm('t','H',H6_new,'P',P6*1000,'CO2');
    S6_new = refpropm('s','T',T6_new,'P',P6*1000,'CO2');
    
    % 更新状态点5和6
    T5 = 0.5*T5 + 0.5*T5_new; % 松弛因子为0.5
    H5 = refpropm('h','T',T5,'P',P5*1000,'CO2');
    S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
    
    T6 = 0.5*T6 + 0.5*T6_new; % 松弛因子为0.5
    H6 = refpropm('h','T',T6,'P',P6*1000,'CO2');
    S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
    
    % 分流点计算
    T12 = T6;
    P12 = P6;
    H12 = H6;
    S12 = S6;
    
    % 副压缩机计算
    S13_is = S12; % 理想等熵过程
    H13_is = refpropm('h','P',P13*1000,'S',S13_is,'CO2');
    H13 = H12 + (H13_is - H12)/N_rc; % 实际焓值
    T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
    S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
    
    % 冷却器计算
    T7 = T6;
    H7 = H6;
    S7 = S6;
    
    % 合流点能量平衡
    H14_calc = alpha*H13 + (1-alpha)*H11;
    
    % 检查收敛
    relative_error = abs((H14_calc - H14_guess)/H14_guess);
    
    if relative_error < 0.001
        object = 0;
    else
        H14_guess = 0.7*H14_guess + 0.3*H14_calc; % 松弛因子为0.3
    end
end

if iteration >= max_iterations
    disp('警告: 达到最大迭代次数，未收敛');
end

% 计算循环效率
% 输入热量
Q_heater = H1 - H16; % 加热器
Q_reheater = H3 - H2; % 再热器
Q_in_total = (Q_heater + Q_reheater)/n_heater;

% 透平输出功率
W_hp_turbine = H1 - H2; % 高压透平
W_lp_turbine = H3 - H4; % 低压透平
W_turbine_total = W_hp_turbine + W_lp_turbine;

% 压缩机消耗功率
W_mc_a = (1-alpha)*(H9 - H8);   % 主压缩机a
W_mc_b = (1-alpha)*(H11 - H10); % 主压缩机b
W_rc = alpha*(H13 - H12);       % 副压缩机
W_comp_total = W_mc_a + W_mc_b + W_rc;

% 净输出功率
W_net = W_turbine_total - W_comp_total;

% 循环热效率
eta_thermal = W_net/Q_in_total;

% 冷却器排出热量
Q_cooler = (1-alpha)*(H7 - H8);
Q_intercooler = (1-alpha)*(H9 - H10);
Q_out_total = Q_cooler + Q_intercooler;

%-----环境温度25℃
T0 = 25 + 273;
H0 = refpropm('h','T',T0,'P',101.325,'CO2');
S0 = refpropm('s','T',T0,'P',101.325,'CO2');

%---火用分析

%---高压透平
E1 = (H1/1000 - H0/1000) - T0*(S1/1000 - S0/1000);
E2 = (H2/1000 - H0/1000) - T0*(S2/1000 - S0/1000);
E_HP_T = (E1 - E2) - (H1/1000 - H2/1000);
N_HP_T = 1-(E_HP_T)/(E1 - E2);

%---再热器
E3 = (H3/1000 - H0/1000) - T0*(S3/1000 - S0/1000);
E_RH = E3 - E2 - (H3/1000 - H2/1000);
N_RH = 1-(E_RH)/(E3 - E2);

%---低压透平
E4 = (H4/1000 - H0/1000) - T0*(S4/1000 - S0/1000);
E_LP_T = (E3 - E4) - (H3/1000 - H4/1000);
N_LP_T = 1-(E_LP_T)/(E3 - E4);

%---高温回热器热侧出口/低温回热器热侧入口 (状态点5)
E5 = (H5/1000 - H0/1000) - T0*(S5/1000 - S0/1000);

%---低温回热器热侧出口/分流点 (状态点6)
E6 = (H6/1000 - H0/1000) - T0*(S6/1000 - S0/1000);

%---主压缩机a
E8 = (H8/1000 - H0/1000) - T0*(S8/1000 - S0/1000);
E9 = (H9/1000 - H0/1000) - T0*(S9/1000 - S0/1000);
E_MC_A = E9 - E8 - (H9/1000 - H8/1000);
N_MC_A = 1-(-E_MC_A)/(E9 - E8);

%---中间冷却器
E10 = (H10/1000 - H0/1000) - T0*(S10/1000 - S0/1000);
E_IC = E9 - E10 - (H9/1000 - H10/1000);
N_IC = 1-(E_IC)/(E9 - E10);

%---主压缩机b
E11 = (H11/1000 - H0/1000) - T0*(S11/1000 - S0/1000);
E_MC_B = E11 - E10 - (H11/1000 - H10/1000);
N_MC_B = 1-(-E_MC_B)/(E11 - E10);

%---副压缩机
E12 = (H12/1000 - H0/1000) - T0*(S12/1000 - S0/1000);
E13 = (H13/1000 - H0/1000) - T0*(S13/1000 - S0/1000);
E_RC = E13 - E12 - (H13/1000 - H12/1000);
N_RC = 1-(-E_RC)/(E13 - E12);

%---冷却器
E7 = (H7/1000 - H0/1000) - T0*(S7/1000 - S0/1000);
E_COOL = E7 - E8 - (H7/1000 - H8/1000);
N_COOL = 1-(E_COOL)/(E7 - E8);

%---低温回热器
E14 = (H14/1000 - H0/1000) - T0*(S14/1000 - S0/1000);
E15 = (H15/1000 - H0/1000) - T0*(S15/1000 - S0/1000);
E_LTR_COLD = E15 - E14 - (H15/1000 - H14/1000); % 冷侧
E_LTR_HOT = E5 - E6 - (H5/1000 - H6/1000);     % 热侧
E_LTR = E_LTR_COLD + E_LTR_HOT;
N_LTR = 1-(E_LTR)/((E15 - E14) + (E5 - E6));

%---高温回热器
E16 = (H16/1000 - H0/1000) - T0*(S16/1000 - S0/1000);
E_HTR_COLD = E16 - E15 - (H16/1000 - H15/1000); % 冷侧
E_HTR_HOT = E4 - E5 - (H4/1000 - H5/1000);     % 热侧
E_HTR = E_HTR_COLD + E_HTR_HOT;
N_HTR = 1-(E_HTR)/((E16 - E15) + (E4 - E5));

%---加热器
E17 = (H17/1000 - H0/1000) - T0*(S17/1000 - S0/1000);
E_HEAT = E17 - E16 - (H17/1000 - H16/1000);
N_HEAT = 1-(E_HEAT)/(E17 - E16);

%---系统总体
E_TOTAL = E_HP_T + E_RH + E_LP_T + E_MC_A + E_IC + E_MC_B + E_RC + E_COOL + E_LTR + E_HTR + E_HEAT;
N_SYSTEM = 1 - E_TOTAL/Q_in_total;

%-----输出结果
disp(['分流比: ',num2str(alpha)])
disp(['THTR: ',num2str(THTR),' K'])
disp(['TLTR: ',num2str(TLTR),' K'])
disp(['循环热效率: ',num2str(eta_thermal*100),' %'])
disp('---火用分析结果---')
disp(['高压透平: ',num2str(E_HP_T),'; 效率=',num2str(N_HP_T*100),'%'])
disp(['再热器: ',num2str(E_RH),'; 效率=',num2str(N_RH*100),'%'])
disp(['低压透平: ',num2str(E_LP_T),'; 效率=',num2str(N_LP_T*100),'%'])
disp(['主压缩机a: ',num2str(-E_MC_A),'; 效率=',num2str(N_MC_A*100),'%'])
disp(['中间冷却器: ',num2str(E_IC),'; 效率=',num2str(N_IC*100),'%'])
disp(['主压缩机b: ',num2str(-E_MC_B),'; 效率=',num2str(N_MC_B*100),'%'])
disp(['副压缩机: ',num2str(-E_RC),'; 效率=',num2str(N_RC*100),'%'])
disp(['冷却器: ',num2str(E_COOL),'; 效率=',num2str(N_COOL*100),'%'])
disp(['低温回热器: ',num2str(E_LTR),'; 效率=',num2str(N_LTR*100),'%'])
disp(['高温回热器: ',num2str(E_HTR),'; 效率=',num2str(N_HTR*100),'%'])
disp(['加热器: ',num2str(E_HEAT),'; 效率=',num2str(N_HEAT*100),'%'])
disp(['系统总体: ',num2str(E_TOTAL),'; 效率=',num2str(N_SYSTEM*100),'%'])