import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pandas.plotting import parallel_coordinates
import seaborn as sns

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来显示中文标签
plt.rcParams['axes.unicode_minus'] = False    # 解决保存图像负号'-'显示为方块的问题

# 读取数据
def load_data(file_path='best.csv'):
    """读取并预处理CSV数据"""
    # 读取数据，没有列名
    df = pd.read_csv(file_path, header=None)
    
    # 根据SCO2_BraytonCycle.m中的信息，为数据添加列名
    columns = [
        'P1', 'P2', 'P4', 'P8', 'alpha', 'n_turb', 'N_mc_a', 'N_mc_b', 
        'neg_eta', 'neg_W_net', 'complexity', 
        'obj4', 'obj5', 'obj6', 'obj7'  # 其他列，可能是约束值或其他信息
    ]
    df.columns = columns
    
    # 转换目标值（因为优化问题中是负号表示最大化）
    df['efficiency'] = -df['neg_eta'] * 100  # 转换为百分比
    df['net_power'] = -df['neg_W_net']      # kJ/kg
    
    return df

# 数据统计分析
def analyze_data(df):
    """基本数据统计分析"""
    print(f"样本数量: {len(df)}")
    
    # 决策变量范围
    print("\n决策变量范围:")
    decision_vars = ['P1', 'P2', 'P4', 'P8', 'alpha', 'n_turb', 'N_mc_a', 'N_mc_b']
    for var in decision_vars:
        print(f"{var}: [{df[var].min():.4f}, {df[var].max():.4f}]")
    
    # 目标函数值范围
    print("\n目标函数值范围:")
    print(f"热效率: [{df['efficiency'].min():.2f}%, {df['efficiency'].max():.2f}%]")
    print(f"净功率: [{df['net_power'].min():.2f} kJ/kg, {df['net_power'].max():.2f} kJ/kg]")
    print(f"复杂度: [{df['complexity'].min():.4f}, {df['complexity'].max():.4f}]")
    
    # Pareto最优解数量
    print(f"\nPareto最优解数量: {len(df)}")
    
    # 典型的Pareto最优解
    print("\n典型Pareto最优解:")
    # 最高效率解
    max_eff_idx = df['efficiency'].idxmax()
    print(f"最高效率解: 效率 = {df.loc[max_eff_idx, 'efficiency']:.2f}%, 净功率 = {df.loc[max_eff_idx, 'net_power']:.2f} kJ/kg, 复杂度 = {df.loc[max_eff_idx, 'complexity']:.4f}")
    
    # 最大功率解
    max_power_idx = df['net_power'].idxmax()
    print(f"最大功率解: 效率 = {df.loc[max_power_idx, 'efficiency']:.2f}%, 净功率 = {df.loc[max_power_idx, 'net_power']:.2f} kJ/kg, 复杂度 = {df.loc[max_power_idx, 'complexity']:.4f}")
    
    # 最小复杂度解
    min_complex_idx = df['complexity'].idxmin()
    print(f"最小复杂度解: 效率 = {df.loc[min_complex_idx, 'efficiency']:.2f}%, 净功率 = {df.loc[min_complex_idx, 'net_power']:.2f} kJ/kg, 复杂度 = {df.loc[min_complex_idx, 'complexity']:.4f}")
    
    # 平衡解
    # 简单方法：找到所有目标函数归一化之和最小的解
    df['norm_eff'] = 1 - (df['efficiency'] - df['efficiency'].min()) / (df['efficiency'].max() - df['efficiency'].min())
    df['norm_power'] = 1 - (df['net_power'] - df['net_power'].min()) / (df['net_power'].max() - df['net_power'].min())
    df['norm_complex'] = (df['complexity'] - df['complexity'].min()) / (df['complexity'].max() - df['complexity'].min())
    df['balance'] = df['norm_eff'] + df['norm_power'] + df['norm_complex']
    balance_idx = df['balance'].idxmin()
    print(f"平衡解: 效率 = {df.loc[balance_idx, 'efficiency']:.2f}%, 净功率 = {df.loc[balance_idx, 'net_power']:.2f} kJ/kg, 复杂度 = {df.loc[balance_idx, 'complexity']:.4f}")
    
    return max_eff_idx, max_power_idx, min_complex_idx, balance_idx

# 绘制决策空间分布
def plot_decision_space(df):
    """绘制决策变量的分布"""
    decision_vars = ['P1', 'P2', 'P4', 'P8', 'alpha', 'n_turb', 'N_mc_a', 'N_mc_b']
    
    plt.figure(figsize=(12, 8))
    plt.boxplot([df[var] for var in decision_vars], labels=decision_vars)
    plt.title('决策变量分布')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig('decision_space_boxplot.png', dpi=300, bbox_inches='tight')
    
    # 散点矩阵
    plt.figure(figsize=(14, 10))
    pd.plotting.scatter_matrix(df[decision_vars], diagonal='kde', alpha=0.7)
    plt.suptitle('决策变量散点矩阵')
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    plt.savefig('decision_space_scatter_matrix.png', dpi=300, bbox_inches='tight')

# 绘制目标空间分布
def plot_objective_space(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx):
    """绘制目标空间的分布"""
    # 3D散点图
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # 所有Pareto解
    ax.scatter(df['efficiency'], df['net_power'], df['complexity'], 
              c='blue', marker='o', alpha=0.6, label='Pareto解')
    
    # 特殊解
    ax.scatter(df.loc[max_eff_idx, 'efficiency'], df.loc[max_eff_idx, 'net_power'], df.loc[max_eff_idx, 'complexity'], 
              c='red', marker='*', s=200, label='最高效率解')
    ax.scatter(df.loc[max_power_idx, 'efficiency'], df.loc[max_power_idx, 'net_power'], df.loc[max_power_idx, 'complexity'], 
              c='green', marker='*', s=200, label='最大功率解')
    ax.scatter(df.loc[min_complex_idx, 'efficiency'], df.loc[min_complex_idx, 'net_power'], df.loc[min_complex_idx, 'complexity'], 
              c='purple', marker='*', s=200, label='最小复杂度解')
    ax.scatter(df.loc[balance_idx, 'efficiency'], df.loc[balance_idx, 'net_power'], df.loc[balance_idx, 'complexity'], 
              c='orange', marker='*', s=200, label='平衡解')
    
    ax.set_xlabel('热效率 (%)')
    ax.set_ylabel('净功率 (kJ/kg)')
    ax.set_zlabel('复杂度')
    ax.set_title('目标空间中的Pareto前沿')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('objective_space_3d.png', dpi=300, bbox_inches='tight')
    
    # 2D投影
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # 效率-功率
    axes[0].scatter(df['efficiency'], df['net_power'], c='blue', alpha=0.6)
    axes[0].scatter(df.loc[max_eff_idx, 'efficiency'], df.loc[max_eff_idx, 'net_power'], c='red', marker='*', s=200)
    axes[0].scatter(df.loc[max_power_idx, 'efficiency'], df.loc[max_power_idx, 'net_power'], c='green', marker='*', s=200)
    axes[0].scatter(df.loc[min_complex_idx, 'efficiency'], df.loc[min_complex_idx, 'net_power'], c='purple', marker='*', s=200)
    axes[0].scatter(df.loc[balance_idx, 'efficiency'], df.loc[balance_idx, 'net_power'], c='orange', marker='*', s=200)
    axes[0].set_xlabel('热效率 (%)')
    axes[0].set_ylabel('净功率 (kJ/kg)')
    axes[0].set_title('效率-功率投影')
    axes[0].grid(True, linestyle='--', alpha=0.7)
    
    # 效率-复杂度
    axes[1].scatter(df['efficiency'], df['complexity'], c='blue', alpha=0.6)
    axes[1].scatter(df.loc[max_eff_idx, 'efficiency'], df.loc[max_eff_idx, 'complexity'], c='red', marker='*', s=200)
    axes[1].scatter(df.loc[max_power_idx, 'efficiency'], df.loc[max_power_idx, 'complexity'], c='green', marker='*', s=200)
    axes[1].scatter(df.loc[min_complex_idx, 'efficiency'], df.loc[min_complex_idx, 'complexity'], c='purple', marker='*', s=200)
    axes[1].scatter(df.loc[balance_idx, 'efficiency'], df.loc[balance_idx, 'complexity'], c='orange', marker='*', s=200)
    axes[1].set_xlabel('热效率 (%)')
    axes[1].set_ylabel('复杂度')
    axes[1].set_title('效率-复杂度投影')
    axes[1].grid(True, linestyle='--', alpha=0.7)
    
    # 功率-复杂度
    axes[2].scatter(df['net_power'], df['complexity'], c='blue', alpha=0.6)
    axes[2].scatter(df.loc[max_eff_idx, 'net_power'], df.loc[max_eff_idx, 'complexity'], c='red', marker='*', s=200)
    axes[2].scatter(df.loc[max_power_idx, 'net_power'], df.loc[max_power_idx, 'complexity'], c='green', marker='*', s=200)
    axes[2].scatter(df.loc[min_complex_idx, 'net_power'], df.loc[min_complex_idx, 'complexity'], c='purple', marker='*', s=200)
    axes[2].scatter(df.loc[balance_idx, 'net_power'], df.loc[balance_idx, 'complexity'], c='orange', marker='*', s=200)
    axes[2].set_xlabel('净功率 (kJ/kg)')
    axes[2].set_ylabel('复杂度')
    axes[2].set_title('功率-复杂度投影')
    axes[2].grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('objective_space_2d_projections.png', dpi=300, bbox_inches='tight')

# 绘制平行坐标图
def plot_parallel_coordinates(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx):
    """绘制平行坐标图"""
    # 创建类别标签
    df['solution_type'] = 'Pareto解'
    df.loc[max_eff_idx, 'solution_type'] = '最高效率解'
    df.loc[max_power_idx, 'solution_type'] = '最大功率解'
    df.loc[min_complex_idx, 'solution_type'] = '最小复杂度解'
    df.loc[balance_idx, 'solution_type'] = '平衡解'
    
    # 归一化数据以便于比较
    cols_to_plot = ['P1', 'P2', 'P4', 'P8', 'alpha', 'n_turb', 'N_mc_a', 'N_mc_b', 
                    'efficiency', 'net_power', 'complexity']
    
    df_norm = df[cols_to_plot].copy()
    for col in cols_to_plot:
        df_norm[col] = (df[col] - df[col].min()) / (df[col].max() - df[col].min())
    
    df_norm['solution_type'] = df['solution_type']
    
    # 绘制平行坐标图
    plt.figure(figsize=(15, 8))
    
    # 先绘制普通Pareto解
    pareto_df = df_norm[df_norm['solution_type'] == 'Pareto解']
    parallel_coordinates(pareto_df, 'solution_type', color='#1f77b4', alpha=0.3)
    
    # 绘制特殊解
    for solution_type, color in [
        ('最高效率解', 'red'),
        ('最大功率解', 'green'),
        ('最小复杂度解', 'purple'),
        ('平衡解', 'orange')
    ]:
        special_df = df_norm[df_norm['solution_type'] == solution_type]
        parallel_coordinates(special_df, 'solution_type', color=color, linewidth=3)
    
    plt.title('决策变量和目标函数的平行坐标图')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('parallel_coordinates.png', dpi=300, bbox_inches='tight')

# 绘制热力图
def plot_correlation_heatmap(df):
    """绘制相关性热力图"""
    cols = ['P1', 'P2', 'P4', 'P8', 'alpha', 'n_turb', 'N_mc_a', 'N_mc_b', 
            'efficiency', 'net_power', 'complexity']
    corr = df[cols].corr()
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5)
    plt.title('决策变量和目标函数之间的相关性')
    plt.tight_layout()
    plt.savefig('correlation_heatmap.png', dpi=300, bbox_inches='tight')

# 为特定解绘制雷达图
def plot_radar_chart(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx):
    """为特殊解绘制雷达图"""
    # 准备数据
    categories = ['热效率', '净功率', '压力比(P1/P4)', '分流比(alpha)', '透平效率', '压缩机效率']
    
    # 计算压力比
    df['pressure_ratio'] = df['P1'] / df['P4']
    
    # 提取数据
    values = {
        '最高效率解': [
            df.loc[max_eff_idx, 'efficiency'] / df['efficiency'].max(),
            df.loc[max_eff_idx, 'net_power'] / df['net_power'].max(),
            df.loc[max_eff_idx, 'pressure_ratio'] / df['pressure_ratio'].max(),
            df.loc[max_eff_idx, 'alpha'] / df['alpha'].max(),
            df.loc[max_eff_idx, 'n_turb'],
            df.loc[max_eff_idx, 'N_mc_a']
        ],
        '最大功率解': [
            df.loc[max_power_idx, 'efficiency'] / df['efficiency'].max(),
            df.loc[max_power_idx, 'net_power'] / df['net_power'].max(),
            df.loc[max_power_idx, 'pressure_ratio'] / df['pressure_ratio'].max(),
            df.loc[max_power_idx, 'alpha'] / df['alpha'].max(),
            df.loc[max_power_idx, 'n_turb'],
            df.loc[max_power_idx, 'N_mc_a']
        ],
        '最小复杂度解': [
            df.loc[min_complex_idx, 'efficiency'] / df['efficiency'].max(),
            df.loc[min_complex_idx, 'net_power'] / df['net_power'].max(),
            df.loc[min_complex_idx, 'pressure_ratio'] / df['pressure_ratio'].max(),
            df.loc[min_complex_idx, 'alpha'] / df['alpha'].max(),
            df.loc[min_complex_idx, 'n_turb'],
            df.loc[min_complex_idx, 'N_mc_a']
        ],
        '平衡解': [
            df.loc[balance_idx, 'efficiency'] / df['efficiency'].max(),
            df.loc[balance_idx, 'net_power'] / df['net_power'].max(),
            df.loc[balance_idx, 'pressure_ratio'] / df['pressure_ratio'].max(),
            df.loc[balance_idx, 'alpha'] / df['alpha'].max(),
            df.loc[balance_idx, 'n_turb'],
            df.loc[balance_idx, 'N_mc_a']
        ]
    }
    
    # 绘制雷达图
    angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
    angles += angles[:1]  # 闭合图形
    
    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    
    for label, val in values.items():
        val_closed = val + val[:1]  # 闭合数据
        ax.plot(angles, val_closed, linewidth=2, label=label)
        ax.fill(angles, val_closed, alpha=0.1)
    
    ax.set_theta_offset(np.pi / 2)  # 旋转图形使第一个轴位于顶部
    ax.set_theta_direction(-1)  # 顺时针方向
    
    # 设置标签
    plt.xticks(angles[:-1], categories)
    
    # 设置雷达图的范围
    ax.set_ylim(0, 1)
    
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.title('特殊解的性能雷达图')
    plt.tight_layout()
    plt.savefig('radar_chart.png', dpi=300, bbox_inches='tight')

# 执行所有分析
def main():
    # 读取数据
    df = load_data()
    
    # 分析数据
    max_eff_idx, max_power_idx, min_complex_idx, balance_idx = analyze_data(df)
    
    # 创建可视化
    plot_decision_space(df)
    plot_objective_space(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx)
    
    try:
        plot_parallel_coordinates(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx)
    except Exception as e:
        print(f"平行坐标图绘制失败: {e}")
    
    plot_correlation_heatmap(df)
    plot_radar_chart(df, max_eff_idx, max_power_idx, min_complex_idx, balance_idx)
    
    plt.close('all')
    print("\n所有图表已保存到当前目录")

if __name__ == "__main__":
    main() 