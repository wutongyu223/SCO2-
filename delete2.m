%%------初始参数
%%---透平
n2 = 0.93;%透平效率
P14 = str2double('20');   %进口压力20，单位Mpa
T14 = str2double('873');  %进口温度，873K;600℃
S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
P1 = str2double('7.6');   %出口压力，单位Mpa
S1 = S14;
H1 = refpropm('h','P',P1*1000,'S',S1,'CO2');
T1 = refpropm('t','H',H1,'P',P1*1000,'CO2');%点1的温度
H1 = H14 + (H1 - H14)*n2; %点1实际的焓
S1 = refpropm('s','P',P1*1000,'H',H1,'CO2');%点1实际的熵

%%---主压缩机（1-x分流，效率N1）
N1 = 0.89;
T5 = str2double('305'); %令主压缩机进口温度32℃
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
P9 = P1;                 %热压缩机进口与透平出口相连，压力相等

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

%X = sdpvar(1,60);
%Y = sdpvar(1,60);
%for x = 1:1:60
x = 0.35;
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
n = ((H14 - H1) - x*(H10 - H9) - (1-x)*(H6 - H5))/((H14 - H12)/n1);
W1 = x*(H10 - H9) + (1-x)*(H6 - H5);
%-----环境温度25℃
T0 = 25 + 273;
H0 = refpropm('h','T',T0,'P',101.325,'CO2');
S0 = refpropm('s','T',T0,'P',101.325,'CO2');
%---透平
E14 = (H14/1000 - H0/1000) - T0*(S14/1000 - S0/1000);
E1 = (H1/1000 - H0/1000) - T0*(S1/1000 - S0/1000);
ET = (E14 - E1) - (H14/1000 - H1/1000);
NT = 1-(ET)/(E14 - E1);
%---主压缩机
E6 = (H6/1000 - H0/1000) - T0*(S6/1000 - S0/1000);
E5 = (H5/1000 - H0/1000) - T0*(S5/1000 - S0/1000);
EZYS = E6 - E5 - (H6/1000 - H5/1000);
NZYS = 1-(-EZYS)/(E6 - E5);
%---副压缩机
E10 = (H10/1000 - H0/1000) - T0*(S10/1000 - S0/1000);
E9 = (H9/1000 - H0/1000) - T0*(S9/1000 - S0/1000);
EFYS = E10 - E9 - (H10/1000 - H9/1000);
NFYS = 1-(-EFYS)/(E10 - E9);
%---冷却器
T4 = T9;
S4 = refpropm('s','T',T4,'H',H4,'CO2');
E4 = (H4/1000 - H0/1000) - T0*(S4/1000 - S0/1000);
E5 = (H5/1000 - H0/1000) - T0*(S5/1000 - S0/1000);
EL = E4 - E5;
NL = E5/E4;
%---高温回热器
S2 = refpropm('s','T',T2,'H',H2,'CO2');
S3 = refpropm('s','T',T3,'H',H3,'CO2');
S12 = refpropm('s','T',T12,'H',H12,'CO2');
S11 = refpropm('s','T',T11,'H',H11,'CO2');
E2 = (H2/1000 - H0/1000) - T0*(S2/1000 - S0/1000);
E3 = (H3/1000 - H0/1000) - T0*(S3/1000 - S0/1000);
E12 = (H12/1000 - H0/1000) - T0*(S12/1000 - S0/1000);
E11 = (H11/1000 - H0/1000) - T0*(S11/1000 - S0/1000);
EHTR = E2 - E3;
NHTR = (E12 - E11)/(E2 - E3);
%---低温回热器
S8 = refpropm('s','T',T8,'H',H8,'CO2');
S7 = refpropm('s','T',T7,'H',H7,'CO2');
E7 = (H7/1000 - H0/1000) - T0*(S7/1000 - S0/1000);
E8 = (H8/1000 - H0/1000) - T0*(S8/1000 - S0/1000);
ELTR = E3 - E4;
NLTR = ((1-x)*(E8 - E7))/(E3 - E4);

disp(['分流比',num2str(x)])
disp(['THTR',num2str(THTR)])
disp(['TLTR',num2str(TLTR)])
disp(['循环效率',num2str(n)])
disp(['透平',num2str(ET),';',num2str(NT)])
disp(['主压缩机',num2str(-EZYS),';',num2str(NZYS)])
disp(['副压缩机',num2str(-EFYS),';',num2str(NFYS)])
disp(['冷却器',num2str(EL),';',num2str(NL)])
disp(['低温回热器',num2str(ELTR),';',num2str(NLTR)])
disp(['高温回热器',num2str(EHTR),';',num2str(NHTR)])