%%------初始参数
%%---透平
n2 = 0.93;%透平效率
P14 = str2double('25');   %进口压力，单位Mpa
T14 = str2double('840');  %进口温度，单位K
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
P1 = str2double('7.5');  %出口压力，单位Mpa
S1 = S14;
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');
H1 = H14 + (H1 - H14)*n2;   %点1实际的焓
T1 = refpropm('t','H',H1,'P',P1*1000,'CO2');
S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');

%%---主压缩机（1-x分流，效率N1）
N1 = 0.89;
T5 = str2double('305'); %令主压缩机进口温度32（效率曲线）
P5 = P1;                %主压缩机进口压力，P5等于透平出口压力P1
S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
H5 = refpropm('h','T',T5,'P',P5*1000,'CO2');
S6 = S5;                %假设压缩完熵不变
P6 = P14;               %冷压缩机出口与透平进口相连，压力相同
H6 = refpropm('h','P',P6*1000,'S',S6,'CO2');
H6 = H5 + (H6 - H5)/N1;
T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');

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
x = 0.32;
for TLTR = 4:1:20
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
        H10 = H9 + (H10 - H9)/N2;%热压缩机实际出口焓值
        T10 = refpropm('t','P',P10*1000,'S',S10,'CO2');
        %冷却器
        H4 = H9;
        %分流汇合
        H11 = x*H10 + (1-x)*H8;
        if (0.98 <= H11/H111) && (H11/H111<= 1.02)
            object = 0;
        else
            H111 = H111 + 50;
        end
    end
    n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1);%循环功：透平-压缩机耗功
    W = 300;%输出净功率
    qm = (W*10^(6))/((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5));%(忽略发电机效率)-计算流量qm，单位kg/s
    W1 = x*(H10 - H9)*qm;
    W2 = (1-x)*(H6 - H5)*qm;
    %高温回热器
    e = exp(1);
    Qh = (H2 - H3)*qm;%单位W
    T212 = T2 - T12;
    T311 = T3 - T11;
    TH = (max(T212,T311) - min(T212,T311))/log(max(T212,T311)/min(T212,T311));
    UH = Qh/TH;
    ftHTR = 1;%最高温度小于550℃，温度修正系数为1
    fpHTR = 1;%压力修正系数
    aHTR = 49.45;%
    bHTR = 0.7544;%
    CHTR = aHTR*((UH)^bHTR)*ftHTR*fpHTR;
    %低温回热器
    T38 = T3 - T8;
    T47 = T9 - T7;
    Ql = (H3 - H9)*qm;%单位W
    TL = (max(T38,T47) - min(T38,T47))/log(max(T38,T47)/min(T38,T47));
    ftLTR = 1;%最高温度小于550℃，温度修正系数为1
    fpLTR = 1;%压力修正系数
    aLTR = 49.45;%
    bLTR = 0.7544;
    UL = Ql/TL;
    CLTR = aLTR*((UL)^bLTR)*ftLTR*fpLTR;

    disp(['TLTR',num2str(TLTR)])
    disp(['循环效率',num2str(n*100)])
    disp(['高温回热器成本',num2str(CHTR*10^(-6))])
    disp(['低温回热器成本',num2str(CLTR*10^(-6))])
end