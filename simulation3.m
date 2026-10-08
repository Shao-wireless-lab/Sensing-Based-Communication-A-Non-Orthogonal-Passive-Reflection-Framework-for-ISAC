%==========================================================================
% Global Optimization for Reverse Engineering System Parameters
% Excludes N=1 from fitting & plotting (Evaluates N = 2 to 10)
% Subplot 1: Per-Tag Data Rate vs N (Mbps)
% Subplot 2: Total Sum Data Rate vs N (Mbps)
% Font Size: 18 throughout
%==========================================================================
clc; clear; close all;

%% 1. Input Target Simulation Data (Ground Truth)
% 第一列为 N=1 的数据 (250 Mbps)
R_sim_M = [
    250, 30.12048, 20.32520, 15.882278, 12.07729, 11.90476, 10.08065,  7.37463, 7.28863, 5.81395;      % 0.5 m
    250, 35.21127, 19.68504, 14.53488,   12.37624, 10.20408,  8.77193,  7.30994, 6.03865, 4.13223;      % 0.7 m
    250, 51.02041, 37.31343, 28.73563,  19.37984, 10.54852,  5.54324,  2.93772, 1.53941, 0.87443       % 1.0 m
];
R_sim = R_sim_M * 1e6; % 转换为 Hz
spacings = [0.5, 0.7, 1.0];
max_N = 10;
N_plot = 2:max_N;      % 从第二个点开始分析与绘图

%% 2. Fixed System Parameters (Known)
B = 250e6 - 10000;      % 有效带宽
Q = 3.0902;             % Q-factor
beta = 0.1;             % 幅度衰减系数
A_n_scale = 0.04;       % 固定 A_n_scale 为 0.04

%% 3. Global Optimization Setup (Multi-Start Search, Excluding N=1)
options = optimset('Display', 'off', 'MaxIter', 3000, 'MaxFunEvals', 3000, 'TolFun', 1e-8);
num_starts = 150;       
best_cost = inf;        
best_x = [];            

fprintf('Running 3D Global Optimization (Excluding N=1, %d starts)...\n', num_starts);
h = waitbar(0, 'Searching for Global Optimum (N=2:10)...');
rng(42); 

for i = 1:num_starts
    waitbar(i/num_starts, h);
    
    % x0 = [xi, P_peak, P_res]
    x0 = [0.5 * rand(), ...        
          1.0 * rand(), ...        
          2.0 * rand()];    
   
    [x_opt, cost] = fminsearch(@(x) bounded_cost_function(x, spacings, max_N, R_sim, B, Q, beta, A_n_scale), x0, options);
    
    if cost < best_cost
        best_cost = cost;
        best_x = x_opt;
    end
end
close(h);

%% 4. Extract Global Optimal Parameters
A_n_opt_scale = A_n_scale;          % 直接使用已知常量 0.04
xi_opt        = abs(best_x(1));     % 从最优解中提取 xi
P_peak_opt    = abs(best_x(2));     % 提取 P_peak
P_res_opt     = abs(best_x(3));     % 提取 P_res
A_n_opt       = A_n_opt_scale / (5^4);

fprintf('\n=================================================\n');
fprintf('*** 3D GLOBAL Optimal Parameters Found (N=2..10) ***\n');
fprintf('A_n scale : %.5e (A_n = %.5e / 5^4) [FIXED]\n', A_n_opt_scale, A_n_opt_scale);
fprintf('xi (SIC)  : %.5f\n', xi_opt);
fprintf('P_peak    : %.5f\n', P_peak_opt);
fprintf('P_res     : %.5f\n', P_res_opt);
fprintf('Minimum Log-MSE Cost: %.4f\n', best_cost);
fprintf('=================================================\n');

%% 5. Generate Theoretical Curves
Rs_theory_opt = zeros(length(spacings), max_N);
for t = 1:length(spacings)
    d_step = spacings(t);
    for N = 1:max_N
        distances = 5 + d_step * (0 : N-1)'; 
        a = 1 ./ (distances.^4); 
        Rs_all = calculate_Rs_original(a, B, Q, beta, P_peak_opt, P_res_opt, A_n_opt, xi_opt, N);
        Rs_theory_opt(t, N) = min(Rs_all);
    end
end

% 提取从第 2 个点开始的数据 (Mbps)
R_sim_plot  = R_sim_M(:, N_plot);
R_th_plot   = Rs_theory_opt(:, N_plot) / 1e6;

% 计算总吞吐量 Total Rate (Mbps) = N * 单标签速率
Total_sim_plot = R_sim_plot .* repmat(N_plot, length(spacings), 1);
Total_th_plot  = R_th_plot .* repmat(N_plot, length(spacings), 1);

%% 6. Approximation Metrics (Goodness of Fit)
% 单标签指标
sim_flat = R_sim_plot(:);
th_flat  = R_th_plot(:);
ss_res = sum((sim_flat - th_flat).^2);
ss_tot = sum((sim_flat - mean(sim_flat)).^2);
R2_single   = 1 - ss_res / ss_tot;
RMSE_single = sqrt(mean((sim_flat - th_flat).^2));
MAPE_single = mean(abs((sim_flat - th_flat) ./ sim_flat)) * 100;

% 总吞吐量指标
tot_sim_flat = Total_sim_plot(:);
tot_th_flat  = Total_th_plot(:);
ss_res_tot = sum((tot_sim_flat - tot_th_flat).^2);
ss_tot_tot = sum((tot_sim_flat - mean(tot_sim_flat)).^2);
R2_total   = 1 - ss_res_tot / ss_tot_tot;
RMSE_total = sqrt(mean((tot_sim_flat - tot_th_flat).^2));
MAPE_total = mean(abs((tot_sim_flat - tot_th_flat) ./ tot_sim_flat)) * 100;

% 控制台打印结果
fprintf('\n--- Global Accuracy Evaluation ---\n');
fprintf('Single Tag Rate: R^2 = %.4f | RMSE = %.3f Mbps | MAPE = %.2f%%\n', R2_single, RMSE_single, MAPE_single);
fprintf('Total Sum Rate : R^2 = %.4f | RMSE = %.3f Mbps | MAPE = %.2f%%\n', R2_total, RMSE_total, MAPE_total);
fprintf('-----------------------------------\n');

%% 7. Visualization: Pure Marker Comparison (No Lines, Font Size = 18)
colors = [
    0.850, 0.325, 0.098;  % 橙色: 0.5 m
    0.466, 0.674, 0.188;  % 绿色: 0.7 m
    0.000, 0.447, 0.741;  % 蓝色: 1.0 m
];
markers = {'s', '^', 'o'};
type_names = {'Spacing = 0.5 m', 'Spacing = 0.7 m', 'Spacing = 1.0 m'};
fig = figure('Name', 'Theoretical vs Simulation Rates', 'Position', [150, 50, 1000, 950], 'Color', 'w');

% --- 子图 1: 单标签速率 (Single Tag Rate) ---
ax1 = subplot(2, 1, 1);
hold on;
% 仿真数据：实心点 (Solid Markers)
for t = 1:length(spacings)
    plot(N_plot, R_sim_plot(t, :), markers{t}, ...
         'LineStyle', 'none', ...
         'Color', colors(t, :), ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', colors(t, :), ...
         'DisplayName', ['Sim: ', type_names{t}]);
end
% 理论数据：同形状/同颜色空心点 (Hollow Markers)
for t = 1:length(spacings)
    plot(N_plot, R_th_plot(t, :), markers{t}, ...
         'LineStyle', 'none', ...
         'Color', colors(t, :), ...
         'LineWidth', 2, ...
         'MarkerSize', 10, ...
         'MarkerFaceColor', 'none', ...
         'DisplayName', ['Theory: ', type_names{t}]);
end
hold off;
xlabel('Number of Tags (N)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('Single Tag Rate (Mbps)', 'FontSize', 18, 'FontName', 'Times New Roman');
title('Per-Tag Data Rate vs N', 'FontSize', 18, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 2, 'FontName', 'Times New Roman');
grid on; set(gca, 'YScale', 'linear', 'FontSize', 18, 'FontName', 'Times New Roman', 'Box', 'on');
xlim([1.5, 10.5]);
xticks(2:10);

% --- 子图 2: 总速率 (Total Sum Rate) ---
ax2 = subplot(2, 1, 2);
hold on;
% 仿真数据：实心点 (Solid Markers)
for t = 1:length(spacings)
    plot(N_plot, Total_sim_plot(t, :), markers{t}, ...
         'LineStyle', 'none', ...
         'Color', colors(t, :), ...
         'MarkerSize', 8, ...
         'MarkerFaceColor', colors(t, :), ...
         'DisplayName', ['Sim: ', type_names{t}]);
end
% 理论数据：同形状/同颜色空心点 (Hollow Markers)
for t = 1:length(spacings)
    plot(N_plot, Total_th_plot(t, :), markers{t}, ...
         'LineStyle', 'none', ...
         'Color', colors(t, :), ...
         'LineWidth', 2, ...
         'MarkerSize', 10, ...
         'MarkerFaceColor', 'none', ...
         'DisplayName', ['Theory: ', type_names{t}]);
end
hold off;
xlabel('Number of Tags (N)', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('Total Rate (Mbps)', 'FontSize', 18, 'FontName', 'Times New Roman');
title('Total System Data Rate vs N', 'FontSize', 18, 'FontName', 'Times New Roman');
legend('Location', 'northeast', 'FontSize', 15, 'NumColumns', 2, 'FontName', 'Times New Roman');
grid on; set(gca, 'YScale', 'linear', 'FontSize', 18, 'FontName', 'Times New Roman', 'Box', 'on');
xlim([1.5, 10.5]);
xticks(2:10);

%% ================== Helper Functions ================== %%
function cost = bounded_cost_function(x, spacings, max_N, R_sim, B, Q, beta, A_n_scale)
    if any(x < 0)
        cost = inf; return;
    end
    
    xi        = x(1);
    P_peak    = x(2);
    P_res     = x(3); 
    A_n       = A_n_scale / (5^4);
    
    cost = 0;
    for t = 1:length(spacings)
        d_step = spacings(t);
        for N = 2:max_N
            distances = 5 + d_step * (0 : N-1)'; 
            a = 1 ./ (distances.^4); 
            
            Rs_all = calculate_Rs_original(a, B, Q, beta, P_peak, P_res, A_n, xi, N);
            Rs_common = min(Rs_all);
            
            log_theory = log10(Rs_common + 1e-10);
            log_sim    = log10(R_sim(t, N) + 1e-10);
            
            cost = cost + (log_theory - log_sim)^2;
        end
    end
end

function Rs_all = calculate_Rs_original(a, B, Q, beta, P_peak, P_res, A_n, xi, N)
    Rs_all = zeros(N, 1);
    for m = 1:N
        sum_prev = 0;
        if m > 1
            sum_prev = xi * sum(a(1:m-1).^2);
        end
        sum_next = 0;
        if m < N
            sum_next = sum(a(m+1:N).^2);
        end
        
        Interf_m = A_n^2 + P_res * ( 0.5 * (1 + beta^2) * (sum_prev + sum_next));
        
        term1 = sqrt(2 * P_peak * a(m)^2 + Interf_m);
        term2 = sqrt(2 * (beta^2) * P_peak * a(m)^2 + Interf_m);
        Denom_m = term1 + term2;
        
        Rs_all(m) = B * ((a(m) * (1 - beta)) / (Q * Denom_m))^2;
    end
end