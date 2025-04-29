classdef SCO2_BraytonCycle < PROBLEM
    % SCO2_BraytonCycle - 超临界CO2布雷顿循环多目标优化问题
    % 该问题可以在PlatEMO中使用NSGA-III算法进行多目标优化
    
    properties
        % 默认参数（用于初始化和范围边界）
        default_P1 = 21;     % 高压透平入口压力(MPa)
        default_P2 = 15;     % 高压透平出口压力(MPa)
        default_P4 = 7.5729; % 低压透平出口压力(MPa)
        default_P8 = 12;     % 主压缩机a出口压力(MPa)
        default_T1 = 873;    % 高压透平入口温度(K)
        default_T7 = 305;    % 冷却器出口温度(K)
        default_alpha = 0.3333; % 默认分流比例
        default_n_hp_turb = 0.93; % 默认高压透平效率
        default_n_lp_turb = 0.93; % 默认低压透平效率
        default_N_mc_a = 0.89; % 默认主压缩机a效率
        default_N_mc_b = 0.89; % 默认主压缩机b效率
        default_N_rc = 0.85;  % 默认副压缩机效率
        default_e1 = 0.86;    % 默认低温回热器效率
        default_e2 = 0.86;    % 默认高温回热器效率
    end
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 3;  % 三个目标：最大化热效率、最大化功率、最小化系统复杂度
            if isempty(obj.D); obj.D = 8; end  % 决策变量数量
            obj.lower = [18, 10, 5, 8, 0.2, 0.85, 0.85, 0.80]; % 决策变量下界
            obj.upper = [25, 18, 10, 15, 0.45, 0.95, 0.95, 0.95]; % 决策变量上界
            obj.encoding = ones(1,obj.D);  % 实数编码
        end

        %% Generate the image of Pareto front
        function R = GetPF(obj)
            R = [];  % 没有已知的前沿
        end
        
        %% 计算目标函数值
        function PopObj = CalObj(obj, PopDec)
            % 获取种群大小
            [N, ~] = size(PopDec);
            % 初始化输出
            PopObj = zeros(N, obj.M);
            
            % 对每个个体计算目标函数值
            for i = 1:N
                % 提取决策变量
                P1 = PopDec(i,1);      % 高压透平入口压力
                P2 = PopDec(i,2);      % 高压透平出口压力
                P4 = PopDec(i,3);      % 低压透平出口压力
                P8 = PopDec(i,4);      % 主压缩机a出口压力
                alpha = PopDec(i,5);   % 分流比例
                n_turb = PopDec(i,6);  % 透平效率
                N_mc = PopDec(i,7);    % 主压缩机效率
                e_recup = PopDec(i,8); % 回热器效率
                
                % 计算循环性能
                try
                    [eta, W_net, complexity] = obj.calcCyclePerformance(...
                        P1, P2, P4, P8, alpha, n_turb, n_turb, N_mc, N_mc, obj.default_N_rc, ...
                        e_recup, e_recup, obj.default_T1, obj.default_T7);
                    
                    % 设置目标：最大化热效率，最大化净功率，最小化复杂度
                    PopObj(i,1) = -eta;                   % 最大化热效率（负号表示最大化）
                    PopObj(i,2) = -W_net/1000;            % 最大化净功率 (kJ/kg)（负号表示最大化）
                    PopObj(i,3) = complexity;             % 最小化系统复杂度
                catch
                    % 如果计算失败，给出一个较差的解
                    PopObj(i,1) = 10;  % 非常低的效率
                    PopObj(i,2) = 10;  % 非常低的净功率
                    PopObj(i,3) = 10;  % 非常高的复杂度
                end
            end
        end
        
        %% 计算约束值
        function PopCon = CalCon(obj, PopDec)
            % 获取种群大小
            [N, ~] = size(PopDec);
            % 初始化约束值（4个约束）
            PopCon = zeros(N, 4);
            
            for i = 1:N
                % 提取决策变量
                P1 = PopDec(i,1);      % 高压透平入口压力
                P2 = PopDec(i,2);      % 高压透平出口压力
                P4 = PopDec(i,3);      % 低压透平出口压力
                P8 = PopDec(i,4);      % 主压缩机a出口压力
                
                % 约束1：确保压力级别合适 P1 > P2 > P4
                PopCon(i,1) = P2 - P1;  % P2 < P1，约束为负
                PopCon(i,2) = P4 - P2;  % P4 < P2，约束为负
                
                % 约束2：确保主压缩机a出口压力小于高压透平入口压力
                PopCon(i,3) = P8 - P1;  % P8 < P1，约束为负
                
                % 约束3：确保主压缩机a出口压力大于低压透平出口压力
                PopCon(i,4) = P4 - P8;  % P4 < P8，约束为负
            end
        end
        
        %% 主要计算函数
        function [eta, W_net, complexity] = calcCyclePerformance(...
                obj, P1, P2, P4, P8, alpha, n_hp_turb, n_lp_turb, N_mc_a, N_mc_b, N_rc, e1, e2, T1, T7)
            % 此函数执行超临界CO2布雷顿循环的核心计算，返回热效率、净功率和复杂度指标
            
            % 固定参数设置
            T3 = T1;               % 再热温度与初始温度相同
            P3 = P2;               % 再热器压力保持不变
            P7 = P4;               % 冷却器出口压力
            P9 = P8;               % 中间冷却压力保持不变
            T9 = T7;               % 中间冷却后温度
            P10 = P1;              % 主压缩机b出口压力
            P11 = P1;              % 副压缩机出口压力
            P6 = P4;               % 分流点压力
            P5 = P4;               % 高温回热器热侧出口压力
            
            %------初始计算------
            % 高压透平入口
            S1 = refpropm('s','T',T1,'P',P1*1000,'CO2');
            H1 = refpropm('h','T',T1,'P',P1*1000,'CO2');
            
            % 高压透平出口
            S2_is = S1;
            H2_is = refpropm('h','P',P2*1000,'S',S2_is,'CO2');
            H2 = H1 - (H1 - H2_is)*n_hp_turb;
            T2 = refpropm('t','P',P2*1000,'H',H2,'CO2');
            S2 = refpropm('s','T',T2,'P',P2*1000,'CO2');
            
            % 再热器出口
            H3 = refpropm('h','T',T3,'P',P3*1000,'CO2');
            S3 = refpropm('s','T',T3,'P',P3*1000,'CO2');
            
            % 低压透平出口
            S4_is = S3;
            H4_is = refpropm('h','P',P4*1000,'S',S4_is,'CO2');
            H4 = H3 - (H3 - H4_is)*n_lp_turb;
            T4 = refpropm('t','P',P4*1000,'H',H4,'CO2');
            S4 = refpropm('s','T',T4,'P',P4*1000,'CO2');
            
            % 冷却器出口/主压缩机a入口
            H7 = refpropm('h','T',T7,'P',P7*1000,'CO2');
            S7 = refpropm('s','T',T7,'P',P7*1000,'CO2');
            
            % 主压缩机a出口
            S8_is = S7;
            H8_is = refpropm('h','P',P8*1000,'S',S8_is,'CO2');
            H8 = H7 + (H8_is - H7)/N_mc_a;
            T8 = refpropm('t','P',P8*1000,'H',H8,'CO2');
            S8 = refpropm('s','T',T8,'P',P8*1000,'CO2');
            
            % 中间冷却器出口
            H9 = refpropm('h','T',T9,'P',P9*1000,'CO2');
            S9 = refpropm('s','T',T9,'P',P9*1000,'CO2');
            
            % 主压缩机b出口
            S10_is = S9;
            H10_is = refpropm('h','P',P10*1000,'S',S10_is,'CO2');
            H10 = H9 + (H10_is - H9)/N_mc_b;
            T10 = refpropm('t','P',P10*1000,'H',H10,'CO2');
            S10 = refpropm('s','T',T10,'P',P10*1000,'CO2');
            
            %------迭代初始化------
            % 使用我们的新方法初始化
            T5_guess = 600;
            
            % 预估副压缩机状态
            T11_est = 320;
            H11_est = refpropm('h','T',T11_est,'P',P6*1000,'CO2');
            S11_est = refpropm('s','T',T11_est,'P',P6*1000,'CO2');
            
            S12_is = S11_est;
            H12_is = refpropm('h','P',P11*1000,'S',S12_is,'CO2');
            H12 = H11_est + (H12_is - H11_est)/N_rc;
            T12 = refpropm('t','P',P11*1000,'H',H12,'CO2');
            S12 = refpropm('s','T',T12,'P',P11*1000,'CO2');
            
            % 假设合流点状态
            H13_ideal = (1-alpha)*H10 + alpha*H12;
            entropy_gen_factor = 0.03;
            entropy_loss = entropy_gen_factor * abs(H12 - H10);
            H13 = H13_ideal - entropy_loss;
            P13 = P10;
            T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
            S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
            
            % 使用效率法估计回热器温度
            P14 = P13;
            T14 = T13 + e1 * (T5_guess - T13);
            H14 = refpropm('h','T',T14,'P',P14*1000,'CO2');
            S14 = refpropm('s','T',T14,'P',P14*1000,'CO2');
            
            % 高温回热器冷侧出口温度
            P15 = P14;
            T15 = T14 + e2 * (T4 - T14);
            H15 = refpropm('h','T',T15,'P',P15*1000,'CO2');
            S15 = refpropm('s','T',T15,'P',P15*1000,'CO2');
            
            % 反向计算高温回热器热侧出口
            Q_HT_cold = H15 - H14;
            H5 = H4 - Q_HT_cold;
            T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
            S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
            
            % 反向计算低温回热器热侧出口
            Q_LT_cold = H14 - H13;
            H6 = H5 - Q_LT_cold;
            T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
            S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
            
            % 更新分流点状态
            T11 = T6;
            H11 = H6;
            S11 = S6;
            
            % 更新副压缩机出口
            S12_is = S11;
            H12_is = refpropm('h','P',P11*1000,'S',S12_is,'CO2');
            H12 = H11 + (H12_is - H11)/N_rc;
            T12 = refpropm('t','P',P11*1000,'H',H12,'CO2');
            S12 = refpropm('s','T',T12,'P',P11*1000,'CO2');
            
            % 重新计算合流点
            H13_ideal = (1-alpha)*H10 + alpha*H12;
            entropy_gen_factor = 0.03;
            entropy_loss = entropy_gen_factor * abs(H12 - H10);
            H13 = H13_ideal - entropy_loss;
            T13 = refpropm('t','P',P13*1000,'H',H13,'CO2');
            S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
            
            %------热平衡迭代求解------
            n_heater = 0.94;
            
            % 热平衡迭代参数
            max_iter = 50;
            relax_coef = 0.5;
            converge_tol = 0.001;
            
            % 迭代求解
            object = 1;
            iter = 0;
            
            while object && iter < max_iter
                iter = iter + 1;
                
                % 保存上一轮值
                H5_old = H5;
                H6_old = H6;
                H12_old = H12;
                H13_old = H13;
                
                % 更新高温回热器
                Q_HT_hot = H4 - H5;
                Q_HT_cold = H13 - H12;
                
                % 调整高温回热器热侧出口焓
                H5_new = H4 - Q_HT_cold;
                
                % 带松弛系数的更新
                relax = relax_coef;
                H5 = H5_old + relax*(H5_new - H5_old);
                T5 = refpropm('t','P',P5*1000,'H',H5,'CO2');
                S5 = refpropm('s','T',T5,'P',P5*1000,'CO2');
                
                % 更新低温回热器
                Q_LT_hot = H5 - H6;
                
                % 计算低温回热器冷侧出口温度
                T10_LTR_out = T10 + e1 * (T5 - T10);
                H10_LTR_out = refpropm('h','T',T10_LTR_out,'P',P10*1000,'CO2');
                S10_LTR_out = refpropm('s','T',T10_LTR_out,'P',P10*1000,'CO2');
                
                % 更新副压缩机入口状态
                T6_rc_in = T6;
                H6_rc_in = H6;
                S6_rc_in = S6;
                
                % 更新副压缩机出口状态
                S11_is = S6_rc_in;
                H11_is = refpropm('h','P',P11*1000,'S',S11_is,'CO2');
                H11 = H6_rc_in + (H11_is - H6_rc_in)/N_rc;
                T11 = refpropm('t','P',P11*1000,'H',H11,'CO2');
                S11 = refpropm('s','T',T11,'P',P11*1000,'CO2');
                
                % 计算低温回热器副路出口温度
                T11_LTR_out = T11 + e1 * (T5 - T11);
                H11_LTR_out = refpropm('h','T',T11_LTR_out,'P',P11*1000,'CO2');
                S11_LTR_out = refpropm('s','T',T11_LTR_out,'P',P11*1000,'CO2');
                
                % 总吸热量
                Q_LT_cold_main = H10_LTR_out - H10;
                Q_LT_cold_bypass = H11_LTR_out - H11;
                Q_LT_cold = (1-alpha)*Q_LT_cold_main + alpha*Q_LT_cold_bypass;
                
                % 调整低温回热器热侧出口焓
                H6_new = H5 - Q_LT_cold;
                H6 = H6_old + relax*(H6_new - H6_old);
                T6 = refpropm('t','P',P6*1000,'H',H6,'CO2');
                S6 = refpropm('s','T',T6,'P',P6*1000,'CO2');
                
                % 计算合流点
                H12_ideal = (1-alpha)*H10_LTR_out + alpha*H11_LTR_out;
                entropy_loss = entropy_gen_factor * abs(H11_LTR_out - H10_LTR_out);
                H12 = H12_ideal - entropy_loss;
                T12 = refpropm('t','P',P10*1000,'H',H12,'CO2');
                S12 = refpropm('s','T',T12,'P',P10*1000,'CO2');
                
                % 更新高温回热器冷侧出口
                T13 = T12 + e2 * (T4 - T12);
                H13 = refpropm('h','T',T13,'P',P13*1000,'CO2');
                S13 = refpropm('s','T',T13,'P',P13*1000,'CO2');
                
                % 检查收敛性
                dH5 = abs(H5 - H5_old)/abs(H5);
                dH6 = abs(H6 - H6_old)/abs(H6);
                dH12 = abs(H12 - H12_old)/abs(H12);
                dH13 = abs(H13 - H13_old)/abs(H13);
                
                if max([dH5, dH6, dH12, dH13]) < converge_tol
                    object = 0;
                end
            end
            
            %------计算循环效率------
            % 功率计算
            W_hp_turb = H1 - H2;
            W_lp_turb = H3 - H4;
            W_mc_a = H8 - H7;
            W_mc_b = H10 - H9;
            W_rc = H11 - H6;
            W_comp_total = (1-alpha)*(W_mc_a + W_mc_b) + alpha*W_rc;
            W_turb_total = W_hp_turb + W_lp_turb;
            W_net = W_turb_total - W_comp_total;
            
            % 热量计算
            Q_heater = H1 - H13;
            Q_reheater = H3 - H2;
            Q_cooler = (1-alpha)*(H6 - H7);
            Q_intercooler = (1-alpha)*(H8 - H9);
            Q_in_total = Q_heater + Q_reheater;
            
            % 循环热效率
            eta = W_net / Q_in_total;
            
            % 复杂度计算 - 定义为一个复合指标
            % 这里使用压力比乘以分流比作为复杂度的简单度量
            pressure_ratio = P1/P4;
            complexity = alpha * pressure_ratio / 10;  % 归一化到合理范围
        end
    end
end 