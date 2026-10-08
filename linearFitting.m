% 清理环境变量
clear; clc; close all;

%% 1. 参数设置
Fs = 50e6;          % 采样率 50 MHz
B = 1e6;            % 信号带宽 1 MHz
% 相关峰的主瓣宽度（用于计算方差时剔除峰值及其附近点）
peak_width = round(Fs / B); 

%% 2. 读取数据
% 读取 CSV 文件中 B22-B1000021 和 C22-C1000021 的数据
% 注意：readmatrix 读取大文件可能需要几秒钟
disp('正在读取CSV文件数据...');
filename = 'tek0001ALL.csv';
data = readmatrix(filename, 'Range', 'B22:C1000021');

ref_sig_all = data(:, 1);  % Reference 信号 (B列)
rx_sig_all = data(:, 2);   % Received 信号 (C列)
disp('数据读取完成！');

% 估算输入信号功率与底噪功率
p_total = var(rx_sig_all); 
% 若信号存在无发射的纯噪声段，可截取该段计算 var(rx_noise)
% 计算理论斜率:
% k_theory = (Fs / B) * (P_sig / P_noise);
% fprintf('理论推导斜率大约应为: %.2f\n', k_theory);

%% 3. 设置变量 T 并进行滑动互相关计算
% 改变 T 等价于改变截取的 Reference 信号长度 (N)
% 这里我们设置一系列不同的 Reference 长度，例如从 500 个点到 10000 个点
ref_lengths = 500:500:100000; 
num_tests = length(ref_lengths);

R0_sq_list = zeros(num_tests, 1);
bg_var_list = zeros(num_tests, 1);
bg_mean_list = zeros(num_tests, 1);

disp('开始进行互相关与统计分析...');
for i = 1:num_tests
    L_ref = ref_lengths(i);
    
    % --- 截取 Reference 信号 ---
    % 从第 1000 个点开始截取 L_ref 个点 (对应你提到的第1000到2000个点等情况)
    ref_start = 1000;
    ref_chunk = ref_sig_all(ref_start : ref_start + L_ref - 1);
    
    % --- 截取 Received 信号 ---
    % 确保 received 包含 reference 并提供足够的滑动范围
    % 设定滑动窗口长度为 2000 个点，这样互相关结果会有 2000+1 个值
    rx_start = 500; 
    L_rx = L_ref + 2000; 
    rx_chunk = rx_sig_all(rx_start : rx_start + L_rx - 1);
    
    % --- 手动互相关 (仅限完全 Overlap 部分) ---
    % 移动 Reference 信号，每次滑动一个点，直到尾部对齐
    num_shifts = L_rx - L_ref + 1;
    R = zeros(num_shifts, 1);
    
    for k = 1:num_shifts
        % 截取 received 中与 reference 长度相等的对应段
        rx_window = rx_chunk(k : k + L_ref - 1);
        % 计算内积作为该位置的互相关值
        R(k) = sum(ref_chunk .* rx_window);
    end
    
    % --- 寻找互相关峰值 R(0) ---
    [peak_val, peak_idx] = max(abs(R)); % 找到绝对值最大的互相关峰
    R0_sq_list(i) = peak_val^2;         % 记录 R(0)^2
    
    % --- 分离背景高斯信号 ---
    % 创建一个掩码，用于剔除峰值所在的主瓣区域
    mask = true(size(R));
    exclude_start = max(1, peak_idx - peak_width);
    exclude_end = min(length(R), peak_idx + peak_width);
    mask(exclude_start:exclude_end) = false;
    
    % 提取除峰值以外的背景序列
    bg_noise = R(mask);
    
    % 计算该 T 下，背景噪声的均值和方差
    bg_mean_list(i) = mean(bg_noise);
    bg_var_list(i) = var(bg_noise);
end
disp('分析完成！');

%% 5. 验证 |R[2l_1]|^2 / \sigma_{bg}^2 与 BT 的关系并计算互相关处理前的输入 SINR
% 1. 计算时间 T 和 时间-带宽积 BT
T_list = (ref_lengths / Fs)';       % 转换为列向量
BT_list = B .* T_list;              % 时间-带宽积 BT (列向量)

% 2. 计算互相关输出指标 (Output SINR)
metric_out = R0_sq_list ./ bg_var_list; % 实际测量值 y

% 3. 强制过零点线性拟合: y = k * x
% 最小二乘求解斜率 k = (x' * y) / (x' * x)
k_fit = BT_list \ metric_out;
fit_y = k_fit * BT_list;            % 拟合值 y_hat（截距严格为 0）

% ==================== 4. 拟合优度计算 ====================
y_real = metric_out;
y_hat  = fit_y;
N = length(y_real);

% 残差平方和 (SSE)
res = y_real - y_hat;
SSE = sum(res.^2);

% 总平方和 (SST_centered: 传统均值中心化; SST_uncentered: 过原点无截距)
SST_centered   = sum((y_real - mean(y_real)).^2);
SST_uncentered = sum(y_real.^2);

% 决定系数 R^2
R2_centered   = 1 - (SSE / SST_centered);     % 常规中心化 R^2
R2_uncentered = 1 - (SSE / SST_uncentered);   % 统计学过原点专用 R^2

% 均方根误差 (RMSE) 与 平均绝对误差 (MAE)
RMSE = sqrt(SSE / (N - 1));                   % 自由度为 N - 1 (估计了1个参数 k)
MAE  = mean(abs(res));

% 提取处理前的原始输入信干噪比
sinr_in_raw_linear = k_fit;
sinr_in_raw_dB = 10 * log10(sinr_in_raw_linear);

% 控制台打印结果与拟合指标
fprintf('\n================== 处理前原始 SINR 评估 (强制过原点拟合) ==================\n');
fprintf('未经互相关处理前的原始输入 SINR (线性值): %.6e\n', sinr_in_raw_linear);
fprintf('未经互相关处理前的原始输入 SINR (dB):     %.2f dB\n', sinr_in_raw_dB);
fprintf('拟合斜率 k:                              %.6e\n', k_fit);
fprintf('拟合截距:                                0.0000 (已强制固定为 0)\n');
fprintf('-------------------------------- 拟合优度指标 --------------------------------\n');
fprintf('常规决定系数 (Centered R^2):             %.4f (%.2f%%)\n', R2_centered, R2_centered * 100);
fprintf('无截距决定系数 (Uncentered R^2):         %.4f (%.2f%%)\n', R2_uncentered, R2_uncentered * 100);
fprintf('均方根误差 (RMSE):                       %.4e\n', RMSE);
fprintf('平均绝对误差 (MAE):                       %.4e\n', MAE);
fprintf('残差平方和 (SSE):                        %.4e\n', SSE);
fprintf('=========================================================================\n\n');

% ---------------------------------------------------------
% 绘图: 符合 IEEE 格式 (图例中带拟合方程及 R^2)
% ---------------------------------------------------------
figure('Color', 'w', 'Units', 'pixels', 'Position', [150, 150, 560, 420]);

% 绘制离散数据点
plot(BT_list, metric_out, 'o', ...
    'LineWidth', 0.5, ...
    'MarkerSize', 2.8, ...
    'MarkerEdgeColor', [0.0, 0.447, 0.741], ...
    'MarkerFaceColor', [0.0, 0.447, 0.741]);
hold on;
plot(BT_list, fit_y, 'r-', 'LineWidth', 1.5);
hold off;

% 坐标轴标签 (严格对应 LaTeX 数学表达)
xlabel('Time-Bandwidth Product $BT$', ...
    'Interpreter', 'latex', 'FontSize', 12);
ylabel('$|R[2l_1]|^2 / \sigma_{\mathrm{bg}}^2$', ...
    'Interpreter', 'latex', 'FontSize', 12);

% 图例说明: 包含 R^2 显示拟合质量
legend_str_fit = sprintf('Linear Fit ($y = kx$, $R^2 = %.4f$)', R2_centered);
legend({'Experimental Measurements', legend_str_fit}, ...
    'Interpreter', 'latex', ...
    'FontSize', 10, ...
    'Location', 'northwest', ...
    'Box', 'on');

% 坐标轴与排版美化
grid on;
set(gca, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'TickLabelInterpreter', 'latex', ...
    'XColor', [0.15, 0.15, 0.15], ...
    'YColor', [0.15, 0.15, 0.15], ...
    'Box', 'on');