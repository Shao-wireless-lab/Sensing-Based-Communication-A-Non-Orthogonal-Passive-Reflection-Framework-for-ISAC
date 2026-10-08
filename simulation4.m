%% 下行通信速率随 Tag 数量增加的演变对比 (BPSK-OFDM 真实物理层)
%clear; clc; close all;

%% 1. 下行 BPSK-OFDM 参数设置
B0 = 250e6;               % 系统总带宽: 250 MHz
B_tag = 3e6;              % 单个 Tag 需求带宽: 3 MHz
K_max = 30;               % Tag 节点最大数量
K_list = 1:1:K_max;

% OFDM 结构参数
N_fft = 1024;             % FFT 点数
ratio_data = 0.85;        % 有效数据子载波比例 (扣除DC、导频与虚子载波)
N_data = round(N_fft * ratio_data);

cp_ratio = 1/8;           % 循环前缀比例 (T_cp / T_u)
bits_per_symbol = 1;      % BPSK 调制: 1 bit/symbol
code_rate = 1;            % 信道编码码率 (无编码设为 1，若有卷积/LDPC可改为 1/2 或 3/4)

% 计算 BPSK-OFDM 下行的基准有效传输速率 (bps)
% R_base = B0 * (N_data/N_fft) * (1 / (1 + cp_ratio)) * bits_per_symbol * code_rate
spec_eff = (N_data / N_fft) * (1 / (1 + cp_ratio)) * bits_per_symbol * code_rate; % bps/Hz
R_dl_nominal = B0 * spec_eff; % 标称下行速率 (~188.89 Mbps)

fprintf('>>> 当前 BPSK-OFDM 下行基准速率为: %.2f Mbps (频谱效率: %.3f bps/Hz)\n', ...
    R_dl_nominal / 1e6, spec_eff);

%% 2. 真实工程与协议开销参数
% --- TDM: 往返传播时延与射频收发切换 ---
c = 3e8;                  % 光速 (m/s)
d_max = 30;               % 最大定位/反射距离: 30 m
T_prop = (2 * d_max) / c; % 往返传播时延: 200 ns
T_trx = 2e-6;             % 射频收发切换与PA建立时间: 2 us
T_guard = T_prop + T_trx; % 单次切换保护间隔 (~2.2 us)

T_frame = 1e-3;           % 调度总帧长: 1 ms
T_tag_data = T_frame * (B_tag / B0); % 单个 Tag 占用数据时隙

% --- FDM: 滤波器滚降保护频带 ---
B_guard = 0.5e6;          % 保护频带: 0.5 MHz / Tag
B_per_tag_fdm = B_tag + B_guard;

%% 3. 三种方案的下行有效速率计算 (Mbps)
% (1) 寄生调制 (Parasitic) - 下行时间/频带零侵占
R_dl_parasitic = (R_dl_nominal * ones(size(K_list))) / 1e6;

% (2) 真实 FDM - 带宽缩减 (假定 OFDM 参数等比缩放)
B_dl_fdm = max(0, B0 - K_list .* B_per_tag_fdm);
R_dl_fdm = (B_dl_fdm .* spec_eff) / 1e6;

% (3) 真实 TDM - 包含 (K+1) 次切换保护间隔的时隙压缩
T_overhead_total = (K_list + 1) .* T_guard;
T_ul_data_total  = K_list .* T_tag_data;
T_dl_available   = max(0, T_frame - T_ul_data_total - T_overhead_total);
R_dl_tdm = ((T_dl_available ./ T_frame) .* R_dl_nominal) / 1e6;

% (4) 理想基准线 (无任何保护时隙与保护频带)
R_dl_ideal = (max(0, 1 - K_list .* (B_tag / B0)) .* R_dl_nominal) / 1e6;

%% 4. 高清绘图输出
figure('Color', [1 1 1], 'Position', [250 200 850 540]);

plot(K_list, R_dl_parasitic, 'r-', 'LineWidth', 2.8, ...
    'DisplayName', 'Parasitic Modulation (Proposed)');
hold on;
plot(K_list, R_dl_ideal, 'k:', 'LineWidth', 1.6, ...
    'DisplayName', 'Ideal TDM / FDM (No Overhead)');
plot(K_list, R_dl_fdm, 'b--s', 'MarkerIndices', 1:5:K_max, 'MarkerSize', 6.5, ...
    'LineWidth', 1.8, 'DisplayName', 'FDM (3 MHz + Guard Band)');
plot(K_list, R_dl_tdm, 'm-.^', 'MarkerIndices', 1:5:K_max, 'MarkerSize', 6.5, ...
    'LineWidth', 1.8, 'DisplayName', 'TDM (3 MHz-equiv. + Guard Interval)');

% 坐标与细节美化
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'GridLineStyle', '--', 'GridAlpha', 0.4);
xlabel('Number of Uplink Tags (K)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Downlink Data Rate (Mbps)', 'FontSize', 14, 'FontWeight', 'bold');
title('Downlink Rate vs. Number of Tags (BPSK-OFDM, B_0 = 250 MHz)', ...
    'FontSize', 15, 'FontWeight', 'bold');

xlim([1, K_max]);
ylim([0, max(R_dl_parasitic) * 1.15]);

legend('Location', 'southwest', 'FontSize', 11, 'Box', 'on', 'Color', [1 1 1 0.9]);

% 动态标注起点速率
text(20, R_dl_parasitic(1) * 0.94, ...
    sprintf('\\leftarrow Parasitic: Constant at %.1f Mbps', R_dl_parasitic(1)), ...
    'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');