function [state, perf] = calculate_cycle(para)
%> @brief 计算超临界CO₂双透平再热再压缩循环的状态点和性能指标
%>
%> @param para struct, 包含以下字段:
%>   - P_high [Pa] 循环最高压力
%>   - P_low [Pa] 循环最低压力
%>   - T_high [K] 循环最高温度
%>   - T_low  [K] 循环最低温度
%>   - P_reheat [Pa] 再热器入口压力
%>   - P_intercool [Pa] 中间冷却器出口压力
%>   - eta_t_HP      等熵效率 (高压透平)
%>   - eta_t_LP      等熵效率 (低压透平)
%>   - eta_c_main    主压缩机等熵效率
%>   - eta_c_recomp  副压缩机等熵效率
%>   - eta_recup_HT  高温回热器效率
%>   - eta_recup_LT  低温回热器效率
%>   - eta_heater    加热器效率
%>   - eta_reheater  再热器热效率
%>   - m_dot [kg/s]  全系统工质质量流量
%>   - alpha         分流比 (副路质量流占比)
%>   - deltaT_HT [K] 高温回热器端差
%>   - deltaT_LT [K] 低温回热器端差
%> @return state struct(17), 每点: T [K], P [Pa], h [kJ/kg], s [kJ/(kg·K)]
%> @return perf  struct, 包含:
%>   - eta_th 热效率
%>   - W_net  净功率 [kW]
%>   - Q_in   热输入 [kW]
%>   - Q_cool 冷却负荷 [kW]
%>   - status 0=OK,1=propertyFail,2=diverge
%> @author 自动生成
%> @date   2024-06-XX

% 默认状态
perf.status = 0;
perf.W_net  = NaN;
perf.Q_in   = NaN;
perf.Q_cool = NaN;
perf.eta_th = NaN;
perf.errorMsg = '';
perf.Q_HT   = NaN;
perf.Q_LT   = NaN;
perf.res_norm = NaN;
% 初始化状态点数组
state = repmat(struct('T',0,'P',0,'h',0,'s',0),17,1);

% 将输入压力从 Pa 转换到 kPa，用于 property 调用
P_high_k      = para.P_high/1e3;
P_low_k       = para.P_low/1e3;
P_reheat_k    = para.P_reheat/1e3;
P_intercool_k = para.P_intercool/1e3;

% 提取其他参数
T_high       = para.T_high;
T_low        = para.T_low;
eta_t_HP     = para.eta_t_HP;
eta_t_LP     = para.eta_t_LP;
eta_c_main   = para.eta_c_main;
eta_c_recomp = para.eta_c_recomp;
eta_recup_HT = para.eta_recup_HT;
eta_recup_LT = para.eta_recup_LT;
eta_heater   = para.eta_heater;
eta_reheater = para.eta_reheater;
m_dot        = para.m_dot;
alpha        = para.alpha;
deltaT_HT    = para.deltaT_HT;
deltaT_LT    = para.deltaT_LT;

max_iter = 80;
tol      = 1e-6;

% 初始猜测：回热器端态
state(16).T = T_high - deltaT_HT; state(16).P = para.P_high;
state(14).T = T_low  + deltaT_LT; state(14).P = para.P_intercool;
% 初始混合点与分流点
state(12).T = state(14).T; state(12).P = state(14).P;
state(7).T  = state(14).T; state(7).P  = para.P_low;
state(15).T = state(14).T; state(15).P = para.P_high;

% 迭代前先初始化所有状态点的焓和熵值
for i = 1:17
    state(i).h = safe_prop('H','T',state(i).T,'P',state(i).P/1e3)/1000;
    state(i).s = safe_prop('S','T',state(i).T,'P',state(i).P/1e3)/1000;
end
% 保存低温回热初始冷侧焓，用于后续计算损失
h14_initial = state(14).h;

try
    for iter = 1:max_iter
        prev_h = [state.h];
        % 1. 高压透平膨胀
        state(1).T = T_high;      state(1).P = para.P_high;
        state(1).h = safe_prop('H','T',state(1).T,'P',state(1).P/1e3)/1000;
        state(1).s = safe_prop('S','T',state(1).T,'P',state(1).P/1e3)/1000;
        % 更新加热器出口状态 state(17) 为 HP 透平入口复制 state(1)
        state(17) = state(1);
        % 状态2: 等熵到 P_reheat
        s2s = state(1).s;
        h2s = safe_prop('H','P',P_reheat_k,'S',s2s*1000)/1000;
        state(2).h = state(1).h - eta_t_HP*(state(1).h - h2s);
        state(2).P = para.P_reheat;
        state(2).T = safe_prop('T','P',state(2).P/1e3,'H',state(2).h*1000);
        state(2).s = safe_prop('S','T',state(2).T,'P',state(2).P/1e3)/1000;
        % 状态3: 再热到最高温度
        state(3).T = T_high; state(3).P = state(2).P;
        state(3).h = safe_prop('H','T',state(3).T,'P',state(3).P/1e3)/1000;
        state(3).s = safe_prop('S','T',state(3).T,'P',state(3).P/1e3)/1000;
        % 低压透平
        s4s = state(3).s;
        h4s = safe_prop('H','P',P_low_k,'S',s4s*1000)/1000;
        state(4).h = state(3).h - eta_t_LP*(state(3).h - h4s);
        state(4).P = para.P_low;
        state(4).T = safe_prop('T','P',state(4).P/1e3,'H',state(4).h*1000);
        state(4).s = safe_prop('S','T',state(4).T,'P',state(4).P/1e3)/1000;
        % 2. 高温回热器显式耦合 (基于效率)
        h4   = state(4).h;
        h15  = state(15).h;
        delta_h_HT = h4 - h15;
        % 更新热侧出口 state(5)
        state(5).h = h4 - eta_recup_HT * delta_h_HT;
        state(5).P = state(4).P;
        state(5).T = safe_prop('T','P',state(5).P/1e3,'H',state(5).h*1000);
        state(5).s = safe_prop('S','T',state(5).T,'P',state(5).P/1e3)/1000;
        % 更新冷侧出口 state(16)
        state(16).h = h15 + eta_recup_HT * delta_h_HT;
        state(16).T = safe_prop('T','P',state(16).P/1e3,'H',state(16).h*1000);
        state(16).s = safe_prop('S','T',state(16).T,'P',state(16).P/1e3)/1000;
        % 3. 低温回热器显式耦合 (基于效率)
        h5        = state(5).h;
        % 与初始冷侧焓计算潜在热差
        delta_h_LT = h5 - h14_initial;
        % 更新热侧出口 state(6)
        state(6).h = h5 - eta_recup_LT * delta_h_LT;
        state(6).P = state(5).P;
        state(6).T = safe_prop('T','P',state(6).P/1e3,'H',state(6).h*1000);
        state(6).s = safe_prop('S','T',state(6).T,'P',state(6).P/1e3)/1000;
        % 更新冷侧出口 state(14)
        state(14).h = h14_initial + eta_recup_LT * delta_h_LT;
        state(14).T = safe_prop('T','P',state(14).P/1e3,'H',state(14).h*1000);
        state(14).s = safe_prop('S','T',state(14).T,'P',state(14).P/1e3)/1000;
        % 4. 分流
        state(7)  = state(6);
        state(12) = state(6);
        % 5. 主路压缩机
        state(8).T = T_low; state(8).P = para.P_low;
        state(8).h = safe_prop('H','T',state(8).T,'P',state(8).P/1e3)/1000;
        state(8).s = safe_prop('S','T',state(8).T,'P',state(8).P/1e3)/1000;
        s9s = state(8).s;
        h9s = safe_prop('H','P',P_intercool_k,'S',s9s*1000)/1000;
        state(9).h = state(8).h + (h9s - state(8).h)/eta_c_main;
        state(9).P = para.P_intercool;
        state(9).T = safe_prop('T','P',state(9).P/1e3,'H',state(9).h*1000);
        state(9).s = safe_prop('S','T',state(9).T,'P',state(9).P/1e3)/1000;
        % 中间冷却
        state(10).T = T_low; state(10).P = state(9).P;
        state(10).h = safe_prop('H','T',state(10).T,'P',state(10).P/1e3)/1000;
        state(10).s = safe_prop('S','T',state(10).T,'P',state(10).P/1e3)/1000;
        % 主压缩机第二级
        s11s = state(10).s;
        h11s = safe_prop('H','P',P_high_k,'S',s11s*1000)/1000;
        state(11).h = state(10).h + (h11s - state(10).h)/eta_c_main;
        state(11).P = para.P_high;
        state(11).T = safe_prop('T','P',state(11).P/1e3,'H',state(11).h*1000);
        state(11).s = safe_prop('S','T',state(11).T,'P',state(11).P/1e3)/1000;
        % 6. 副路压缩机
        state(13) = state(12);
        s13s = state(13).s;
        h13s = safe_prop('H','P',P_high_k,'S',s13s*1000)/1000;
        state(13).h = state(12).h + (h13s - state(12).h)/eta_c_recomp;
        state(13).P = para.P_high;
        state(13).T = safe_prop('T','P',state(13).P/1e3,'H',state(13).h*1000);
        state(13).s = safe_prop('S','T',state(13).T,'P',state(13).P/1e3)/1000;
        % 7. 合流点
        state(15).h = (1-alpha)*state(11).h + alpha*state(13).h;
        state(15).P = para.P_high;
        state(15).T = safe_prop('T','P',state(15).P/1e3,'H',state(15).h*1000);
        state(15).s = safe_prop('S','T',state(15).T,'P',state(15).P/1e3)/1000;
        % 收敛检查
        if max(abs([state.h] - prev_h)) < tol
            break;
        end
    end
    if iter == max_iter && max(abs([state.h] - prev_h)) >= tol
        perf.status = 2;
    end
    % 性能指标计算
    W_turb = m_dot*((state(1).h - state(2).h) + (state(3).h - state(4).h));
    % 按质量流分配计算压缩机功
    m_main = (1-alpha) * m_dot;
    m_recomp = alpha * m_dot;
    W_comp = m_main * ((state(9).h - state(8).h) + (state(11).h - state(10).h)) + m_recomp * (state(13).h - state(12).h);
    W_net  = W_turb - W_comp;
    % 计算热输入 Q_in：主加热器（除以效率）和再热器
    Q_heater        = m_dot * (state(17).h - state(16).h);
    Q_reheat        = m_dot * (state(3).h - state(2).h);
    Q_fuel          = Q_heater / eta_heater + Q_reheat / eta_reheater;
    % 回热器不可逆损失（诊断用）
    Q_loss_HT       = (1 - eta_recup_HT) * m_dot * (state(4).h - state(15).h);
    % 低温回热未回收热：基于初始冷侧焓计算损失
    Q_loss_LT       = (1 - eta_recup_LT) * m_dot * (state(5).h - h14_initial);
    % 中间冷却器冷却负荷：仅主路流量
    Q_intercool     = m_main * (state(7).h - state(8).h);
    % 冷却负荷 Q_cool 包含中间冷却与回热损失
    Q_cool          = Q_intercool + Q_loss_HT + Q_loss_LT;
    % 填充 perf 字段
    perf.W_net      = W_net;
    perf.Q_in       = Q_fuel;
    perf.Q_cool     = Q_cool;
    perf.eta_th     = perf.W_net / perf.Q_in;
    perf.Q_loss_HT  = Q_loss_HT;
    perf.Q_loss_LT  = Q_loss_LT;
    perf.res_norm   = max(abs([state.h] - prev_h));
    % 打印能量平衡误差
    perf.eb_error   = abs((perf.Q_in - (perf.W_net + perf.Q_cool)) / perf.Q_in);
    fprintf('Corrected energy error = %.3f%%\n', perf.eb_error*100);
    fprintf('LT Debug: h5_end=%.3f, h14_initial=%.3f, delta_h0=%.3f, Q_loss_LT=%.3f\n', state(5).h, h14_initial, state(5).h - h14_initial, Q_loss_LT);
    fprintf('Detail: Q_heater=%.3f, Q_reheat=%.3f, Q_loss_HT=%.3f, Q_loss_LT=%.3f, Q_intercool=%.3f\n', Q_heater, Q_reheat, Q_loss_HT, Q_loss_LT, Q_intercool);
catch ME
    perf.errorMsg = ME.message;
    if contains(ME.message, 'FLSH')
        perf.status = 1;
    else
        perf.status = 2;
    end
end

% 内部安全的 property 调用
function val = safe_prop(prop, spec1, v1, spec2, v2)
    try
        val = refpropm(prop, spec1, v1, spec2, v2, 'CO2');
    catch err
        if contains(err.message, 'FLSH')
            dv = 1e-6 * max(abs(v1),1);
            v1p = refpropm(prop, spec1, v1 - dv, spec2, v2, 'CO2');
            v2p = refpropm(prop, spec1, v1 + dv, spec2, v2, 'CO2');
            val = (v1p + v2p) / 2;
        else
            rethrow(err);
        end
    end
end

%{
% 单元测试示例:
% para.P_high      = 16e6; para.P_low       = 7.5e6;
% para.T_high      = 840;   para.T_low       = 305;
% para.P_reheat    = 10e6;  para.P_intercool = 10e6;
% para.eta_t_HP    = 0.93;  para.eta_t_LP    = 0.93;
% para.eta_c_main  = 0.89;  para.eta_c_recomp= 0.89;
% para.eta_recup_HT= 0.86;  para.eta_recup_LT= 0.86;
% para.eta_heater  = 0.94;  para.m_dot       = 100;
% para.alpha       = 0.3;   para.deltaT_HT   = 10;
% para.deltaT_LT   = 10;
% REF.T0 = 298.15; REF.P0 = 101.325;
% REF.h0 = refpropm('H','T',REF.T0,'P',REF.P0,'CO2');
% REF.s0 = refpropm('S','T',REF.T0,'P',REF.P0,'CO2');
% [state, perf] = calculate_cycle(para, REF);
% disp(perf);
%}
end 