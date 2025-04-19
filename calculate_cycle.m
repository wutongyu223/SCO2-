function [state,perf] = calculate_cycle(para)
    % CALCULATE_CYCLE  Steady‑state model of a reheated–recompressed
    %                  double‑turbine sCO2 Brayton cycle.
    %
    % ─────────────── 输入参数 ───────────────
    %   para : struct  (新增 deltaT_inter)
    %       P_high  [Pa]  – Max cycle pressure
    %       P_low   [Pa]  – Min cycle pressure
    %       T_high  [K]   – Turbine‑inlet temperature
    %       T_low   [K]   – Compressor‑suction temperature
    %       P_reheat      – Reheater inlet pressure
    %       P_intercool   – Intercooler outlet pressure (= main‑C1 discharge)
    %       deltaT_inter  – Approach temperature of intercooler  (K)
    %       eta_t_HP, eta_t_LP
    %       eta_c_main,  eta_c_recomp
    %       eta_recup_HT, eta_recup_LT
    %       eta_heater,   eta_reheater
    %       m_dot   [kg/s] – Total mass flow
    %       alpha          – Recompression mass fraction
    %       deltaT_HT, deltaT_LT
    %
    % ─────────────── 输出 ───────────────
    %   state(17) : struct array, fields {T,P,h,s}
    %   perf      : struct –   W_net, Q_in, Q_cool, eta_th, Eb_error, status
    %
    % 流程顺序
    %   1  HP‑Turb → 2 → Reheat → 3 → LP‑Turb → 4
    %   4 → HT‑Recup(hot) → 5 → LT‑Recup(hot) → 6 → split α
    %   main‑path : 6 → 7(C1‑in) → 8(C1‑out) → Intercooler → 9 → 10(C2‑out)
    %   recomp    : 6 → 11(RC‑in) → 12(RC‑out)
    %   Merge → 13 → LT‑Recup(cold) → 14 → HT‑Recup(cold) → 15 → Heater → 16(=1)
    %
    %   All heat loads Q = m·Δh,   Q = max(0,Q)  (non‑negative by definition)
    %   Recuperator losses referenced to the *initial cold‑side enthalpy* (h14_0).
    %
    % Author : <your name>      2025‑04‑18
    % ---------------------------------------------------------------
    
    %% ───── 1. 前置工具 ─────
    safeH = @(varargin) safe_prop('H',varargin{:})/1000;   % kJ/kg
    safeS = @(varargin) safe_prop('S',varargin{:})/1000;   % kJ/kg‐K
    safeT = @(varargin) safe_prop('T',varargin{:});
    KPa   = @(P) P/1e3;
    
    nz = @(x) max(0,x);          % non‑negative helper
    
    %% ───── 2. 状态数组预分配 ─────
    state = repmat(struct('T',0,'P',0,'h',0,'s',0),17,1);
    
    % 为可读性定义索引常量
    HP_IN=1;  HP_OUT=2;  RH_OUT=3;  LP_OUT=4;   ...
    HT_hot_out=5;  LT_hot_out=6;  C1_IN=7;  C1_OUT=8;
    IC_OUT=9;       C2_OUT=10; RC_IN=11;  RC_OUT=12;
    MIX=13; LT_cold_out=14; HT_cold_out=15;
    HEATER_IN=16; HEATER_OUT=17;
    
    %% ───── 3. 基本参数与便利变量 ─────
    P_hi = para.P_high;  P_lo = para.P_low;
    P_re = para.P_reheat; P_ic = para.P_intercool;
    T_hi = para.T_high;  T_lo = para.T_low;
    
    mdot = para.m_dot;  a = para.alpha;
    
    % 关键效率
    eta_HP = para.eta_t_HP;  eta_LP = para.eta_t_LP;
    eta_C  = para.eta_c_main; eta_RC = para.eta_c_recomp;
    eps_HT = para.eta_recup_HT; eps_LT = para.eta_recup_LT;
    
    %% ───── 4. ―― Turbines & Reheat ―― ─────
    state(HP_IN).T = T_hi; state(HP_IN).P = P_hi;
    state(HP_IN).h = safeH('T',T_hi,'P',KPa(P_hi));
    state(HP_IN).s = safeS('T',T_hi,'P',KPa(P_hi));
    
    % 等熵出口
    h2s = safeH('P',KPa(P_re),'S',state(HP_IN).s*1000);
    state(HP_OUT).h = state(HP_IN).h - eta_HP*(state(HP_IN).h - h2s);
    state(HP_OUT).P = P_re;
    state(HP_OUT).T = safeT('P',KPa(P_re),'H',state(HP_OUT).h*1000);
    state(HP_OUT).s = safeS('T',state(HP_OUT).T,'P',KPa(P_re));
    
    % 再热 → RH_OUT
    state(RH_OUT).T = T_hi;
    state(RH_OUT).P = state(HP_OUT).P;
    state(RH_OUT).h = safeH('T',T_hi,'P',KPa(P_re));
    state(RH_OUT).s = safeS('T',T_hi,'P',KPa(P_re));
    
    % LP Turbine
    h4s = safeH('P',KPa(P_lo),'S',state(RH_OUT).s*1000);
    state(LP_OUT).h = state(RH_OUT).h - eta_LP*(state(RH_OUT).h - h4s);
    state(LP_OUT).P = P_lo;
    state(LP_OUT).T = safeT('P',KPa(P_lo),'H',state(LP_OUT).h*1000);
    state(LP_OUT).s = safeS('T',state(LP_OUT).T,'P',KPa(P_lo));
    
    %% ───── 5. ―― HT Recuperator ―― ─────
    % 冷侧初始焓在合流前尚未知，用占位符延后迭代
    state(MIX).h = NaN;  % will be overwritten later
    
    % 初始估计：假设 HT effectiveness 达到设计值一次收敛
    d_h_HTmax = state(LP_OUT).h - state(MIX).h;    % MIX.h 待会儿赋值
    % 用匿名函数在合流后再更新
    
    %% ───── 6. ―― LT Recuperator & Split ――
    % 先设一个 LT 冷侧初始温度 = T_lo + deltaT_LT
    state(LT_cold_out).T = T_lo + para.deltaT_LT;
    state(LT_cold_out).P = P_ic;                     % 低温侧初始压力 = intercool出口
    state(LT_cold_out).h = safeH('T',state(LT_cold_out).T,'P',KPa(P_ic));
    state(LT_cold_out).s = safeS('T',state(LT_cold_out).T,'P',KPa(P_ic));
    
    h14_0 = state(LT_cold_out).h;    % *固定* 作为损失基准
    
    %% ───── 7. ―― Main compressor stage‑1 ――
    state(C1_IN) = state(LT_hot_out);   % 位置保留，稍后更新
    % 为防止循环依赖，先把 hot‑LT 出口临时设成 LP_OUT
    state(LT_hot_out) = state(LP_OUT);
    
    % 迭代上限与收敛容差
    for it = 1:60
        prev = [state.h];       %#ok<AGROW>
    
        % 根据最新 LT_hot_out 更新 main‑path入口
        state(C1_IN) = state(LT_hot_out);
    
        % ── C1 (η_c_main)
        s7 = state(C1_IN).s;
        h8s = safeH('P',KPa(P_ic),'S',s7*1000);
        state(C1_OUT).h = state(C1_IN).h + (h8s-state(C1_IN).h)/eta_C;
        state(C1_OUT).P = P_ic;
        state(C1_OUT).T = safeT('P',KPa(P_ic),'H',state(C1_OUT).h*1000);
        state(C1_OUT).s = safeS('T',state(C1_OUT).T,'P',KPa(P_ic));
    
        % ── Intercooler (approach ΔT = deltaT_inter)
        state(IC_OUT).T = T_lo + para.deltaT_inter;
        state(IC_OUT).P = state(C1_OUT).P;
        state(IC_OUT).h = safeH('T',state(IC_OUT).T,'P',KPa(P_ic));
        state(IC_OUT).s = safeS('T',state(IC_OUT).T,'P',KPa(P_ic));
    
        % ── C2
        s9 = state(IC_OUT).s;
        h10s = safeH('P',KPa(P_hi),'S',s9*1000);
        state(C2_OUT).h = state(IC_OUT).h + (h10s-state(IC_OUT).h)/eta_C;
        state(C2_OUT).P = P_hi;
        state(C2_OUT).T = safeT('P',KPa(P_hi),'H',state(C2_OUT).h*1000);
        state(C2_OUT).s = safeS('T',state(C2_OUT).T,'P',KPa(P_hi));
    
        %% ───── 8. ―― Recompressor path ――
        state(RC_IN) = state(LT_hot_out);  % same inlet
        s11 = state(RC_IN).s;
        h12s = safeH('P',KPa(P_hi),'S',s11*1000);
        state(RC_OUT).h = state(RC_IN).h + (h12s-state(RC_IN).h)/eta_RC;
        state(RC_OUT).P = P_hi;
        state(RC_OUT).T = safeT('P',KPa(P_hi),'H',state(RC_OUT).h*1000);
        state(RC_OUT).s = safeS('T',state(RC_OUT).T,'P',KPa(P_hi));
    
        %% ───── 9. ―― Merge point ――
        state(MIX).h = (1-a)*state(C2_OUT).h + a*state(RC_OUT).h;
        state(MIX).P = P_hi;
        state(MIX).T = safeT('P',KPa(P_hi),'H',state(MIX).h*1000);
        state(MIX).s = safeS('T',state(MIX).T,'P',KPa(P_hi));
    
        %% ─────10. ―― Recuperators (now with MIX updated) ――
        % ------- HT -------
        d_h_HTmax = state(LP_OUT).h - state(MIX).h;
        d_h_HT     = eps_HT * d_h_HTmax;
        state(HT_hot_out).h  = state(LP_OUT).h - d_h_HT;
        state(HT_cold_out).h = state(MIX).h   + d_h_HT;
    
        % property calls
        state(HT_hot_out).P = P_lo;
        state(HT_hot_out).T = safeT('P',KPa(P_lo),'H',state(HT_hot_out).h*1000);
        state(HT_hot_out).s = safeS('T',state(HT_hot_out).T,'P',KPa(P_lo));
    
        state(HT_cold_out).P = P_hi;
        state(HT_cold_out).T = safeT('P',KPa(P_hi),'H',state(HT_cold_out).h*1000);
        state(HT_cold_out).s = safeS('T',state(HT_cold_out).T,'P',KPa(P_hi));
    
        % ------- LT -------
        d_h_LTmax = state(HT_hot_out).h - h14_0;
        d_h_LT    = eps_LT * d_h_LTmax;
        state(LT_hot_out).h = state(HT_hot_out).h - d_h_LT;
        state(LT_cold_out).h = h14_0 + d_h_LT;  % cold‑side outlet
    
        state(LT_hot_out).P = P_lo;
        state(LT_hot_out).T = safeT('P',KPa(P_lo),'H',state(LT_hot_out).h*1000);
        state(LT_hot_out).s = safeS('T',state(LT_hot_out).T,'P',KPa(P_lo));
    
        state(LT_cold_out).T = safeT('P',KPa(P_ic),'H',state(LT_cold_out).h*1000);
        state(LT_cold_out).s = safeS('T',state(LT_cold_out).T,'P',KPa(P_ic));
    
        %% ---- 收敛检查 ----
        if max(abs([state.h]-prev)) < 1e-6,    break,   end
    end
    
    %% ───── 11. Heater ─────
    state(HEATER_IN) = state(HT_cold_out);
    state(HEATER_OUT) = state(HP_IN);  % same as turbine inlet
    
    %% ───── 12. 性能计算 ─────
    m_main = (1-a)*mdot;   m_rc = a*mdot;
    
    W_turb = mdot * ( state(HP_IN).h - state(HP_OUT).h  ...
                    + state(RH_OUT).h - state(LP_OUT).h );
    
    W_comp = m_main*( state(C1_OUT).h - state(C1_IN).h ...
                    + state(C2_OUT).h - state(IC_OUT).h ) ...
           + m_rc  *( state(RC_OUT).h - state(RC_IN).h );
    
    perf.W_net = W_turb - W_comp;
    
    % 热输入（加热器+再热器，含效率）
    Q_heater = nz( mdot*(state(HEATER_OUT).h - state(HEATER_IN).h) );
    Q_reheat = nz( mdot*(state(RH_OUT).h     - state(HP_OUT).h   ) );
    perf.Q_in = Q_heater/para.eta_heater + Q_reheat/para.eta_reheater;
    
    % 冷却负荷：IC + 未回收损失
    Q_IC  = nz( m_main*(state(C1_OUT).h - state(IC_OUT).h) );
    Q_LTh = nz( mdot*(d_h_LTmax - d_h_LT) );   % LT 未回收
    Q_HTh = nz( mdot*(d_h_HTmax - d_h_HT) );   % HT 未回收
    perf.Q_cool = Q_IC + Q_LTh + Q_HTh;
    
    perf.eta_th = perf.W_net / perf.Q_in;
    perf.Eb_error = abs( perf.Q_in - (perf.W_net + perf.Q_cool) ) / perf.Q_in;
    perf.status = 0;
    
    %% ───── 13. 打印/返回 ─────
    fprintf('[Self‑test]  Energy error = %.3f %%\n',perf.Eb_error*100);
    
    %% ───── 14. 内置单元测试 (仅首次调用演示) ─────
    if nargin==0  % 若直接运行文件 → demo
        para = default_para();          %#ok<NASGU>
        [~,~] = calculate_cycle(para);  % 递归调用演示
    end
    
    % ====================================================================== %
    function val = safe_prop(prop,varargin)
    % 安全封装 REPROP 调用，若闪蒸失败则尝试微调求平均
        try
            val = refpropm(prop,varargin{:},'CO2');
        catch
            x = varargin{2};  dv = 1e-4*max(abs(x),1);
            v1 = refpropm(prop,varargin{:},'CO2','T',varargin{2}-dv);
            v2 = refpropm(prop,varargin{:},'CO2','T',varargin{2}+dv);
            val = (v1+v2)/2;
        end
    end
    % ---------------------------------------------------------------------- %
    function p = default_para()
    % 默认测试参数（示例）
        p.P_high = 16e6;   p.P_low = 7.5e6;
        p.T_high = 840;    p.T_low = 305;
        p.P_reheat = 10e6; p.P_intercool = 10e6;
        p.deltaT_inter = 8;      % 新增 intercool 端差
        p.eta_t_HP = .93;  p.eta_t_LP = .93;
        p.eta_c_main = .89; p.eta_c_recomp = .89;
        p.eta_recup_HT = .86; p.eta_recup_LT = .86;
        p.eta_heater = .94; p.eta_reheater = .94;
        p.m_dot = 100;     p.alpha = 0.3;
        p.deltaT_HT = 10;  p.deltaT_LT = 10;
    end
    
    % 如果没有输入参数，则使用默认参数
    if nargin < 1
        para = default_para();
    end
    