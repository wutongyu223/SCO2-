%% RunNSGAIII_SCO2.m - 运行NSGA-III算法优化超临界CO2布雷顿循环
% 此脚本用于直接从MATLAB命令行运行NSGA-III算法优化SCO2布雷顿循环
% 需要将SCO2_BraytonCycle.m放在正确的位置

% 检查PlatEMO路径
if ~exist('platemo', 'file')
    error('未找到PlatEMO平台，请确保已安装并添加到MATLAB路径中');
end

% 检查REFPROP是否可用
try
    refpropm('h','T',300,'P',10000,'CO2');
    disp('REFPROP可用，继续执行...');
catch
    error('REFPROP不可用，请确保已安装并配置正确');
end

%% 设置算法参数
algorithm   = 'NSGAIII';     % 算法: NSGA-III
problem     = 'SCO2_BraytonCycle'; % 问题: 超临界CO2布雷顿循环
M           = 3;             % 目标数: 3
popsize     = 91;            % 种群大小: 91 (建议为 (M+H-1)!/(H!*(M-1)!), H=4)
maxFE       = 10000;         % 最大评价次数
encoding    = 'real';        % 编码方式: 实数编码

%% 其他参数
save_interval = 100;         % 保存间隔（每evaluations次评估）
output_path = './Results/';  % 结果输出路径
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

%% 运行算法
disp('开始运行NSGA-III算法优化超临界CO2布雷顿循环...');
tic;
[Population, Obj, ~] = platemo('algorithm', algorithm, ...
                              'problem', problem, ...
                              'M', M, ...
                              'maxFE', maxFE, ...
                              'save', save_interval, ...
                              'outputFcn', @(Algorithm) save_results(Algorithm, output_path), ...
                              'N', popsize);
runtime = toc;

%% 输出结果
fprintf('优化完成，用时 %.2f 秒\n', runtime);
fprintf('找到 %d 个非支配解\n', length(Population));

% 保存最终结果
save([output_path, 'final_results.mat'], 'Population', 'Obj', 'runtime');
saveas(gcf, [output_path, 'pareto_front.fig']);
saveas(gcf, [output_path, 'pareto_front.png']);

% 绘制Pareto前沿
figure('Position', [100, 100, 1200, 400]);

% 第一个子图: 热效率 vs 净功率
subplot(1, 3, 1);
scatter(-Obj(:,1), -Obj(:,2), 30, Obj(:,3), 'filled');
xlabel('热效率'); ylabel('净功率 (kJ/kg)');
title('热效率 vs 净功率');
colorbar; colormap('jet'); 
grid on;

% 第二个子图: 热效率 vs 复杂度
subplot(1, 3, 2);
scatter(-Obj(:,1), Obj(:,3), 30, -Obj(:,2), 'filled');
xlabel('热效率'); ylabel('系统复杂度');
title('热效率 vs 系统复杂度');
colorbar; colormap('jet');
grid on;

% 第三个子图: 净功率 vs 复杂度
subplot(1, 3, 3);
scatter(-Obj(:,2), Obj(:,3), 30, -Obj(:,1), 'filled');
xlabel('净功率 (kJ/kg)'); ylabel('系统复杂度');
title('净功率 vs 系统复杂度');
colorbar; colormap('jet');
grid on;

% 调整并保存图形
sgtitle('超临界CO2布雷顿循环多目标优化结果');
saveas(gcf, [output_path, 'objectives_scatter.fig']);
saveas(gcf, [output_path, 'objectives_scatter.png']);

% 输出最优解详情
fprintf('\n最优解详情：\n');
fprintf('==================================\n');
fprintf('ID\t热效率\t\t净功率(kJ/kg)\t系统复杂度\n');
for i = 1:min(10, length(Population))
    fprintf('%d\t%.4f\t\t%.2f\t\t%.4f\n', i, -Obj(i,1), -Obj(i,2), Obj(i,3));
end

% 显示最高效率和最大功率的解
[~, idx_max_eff] = min(Obj(:,1));
[~, idx_max_power] = min(Obj(:,2));
fprintf('\n最高效率解：效率 = %.4f, 功率 = %.2f kJ/kg, 复杂度 = %.4f\n', ...
        -Obj(idx_max_eff,1), -Obj(idx_max_eff,2), Obj(idx_max_eff,3));
fprintf('最大功率解：效率 = %.4f, 功率 = %.2f kJ/kg, 复杂度 = %.4f\n', ...
        -Obj(idx_max_power,1), -Obj(idx_max_power,2), Obj(idx_max_power,3));

fprintf('\n结果已保存至 %s\n', output_path);

% 显示一些建议
disp('建议：');
disp('1. 若优化结果不理想，可增加评价次数和种群大小');
disp('2. 可以调整决策变量范围以探索更广阔的解空间');
disp('3. 对于特定应用场景，可根据需求选择合适权衡点');

% 结果保存函数
function save_results(Algorithm, output_path)
    if mod(Algorithm.FE, 1000) == 0
        fprintf('已完成评价次数：%d / %d\n', Algorithm.FE, Algorithm.maxFE);
        
        % 保存当前结果
        Population = Algorithm.result{1};
        Obj = Population.objs;
        gen = floor(Algorithm.FE / Algorithm.N);
        
        save([output_path, 'gen_', num2str(gen), '.mat'], 'Population', 'Obj');
    end
end 