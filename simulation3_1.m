% Comprehensive Simulation: Single-Tag Bit Rate & Total Capacity vs. Tag Count (N = 1 to N_max)
% Using Log-Space Binary Search for High Precision SPB Determination
% Comparison across three spacing configurations: 0.5m, 0.7m, and 1.0m
% Baseline Noise is strictly defined by the 5m single-tag echo power

clc; close all; clear;

%% --- 1. System & Physical Channel Configurations ---
N_max = 10;                   % 最大 Tag 数量
target_BER = 0.001;           % 目标 BER 上限 (0.1%)
snr = 15;                     % 基于 5m 单 Tag 回波计算的固定信噪比 (dB)
c = 3e8;                      % 电磁波传播速率 (m/s)

% 间距配置
spacing_list = [0.5, 0.7, 1.0]; % 三组对比间距 (m)
num_spacings = length(spacing_list);

% 信号物理参数
fs = 2.5e9;                   % 采样率 2.5 GHz
t_duration = 0.004;           % 信号持续时间 (s)
f_start = 10000;              % 10 kHz
f_end = 250e6;                % 250 MHz
amp_on = 1;                   % OOK '1' 幅度系数
amp_off = 0.1;                % OOK '0' 幅度系数

%% --- 2. Transmitted OFDM Base Waveform Generation ---
N_fft = 5120;
df = fs / N_fft;
cp_ratio = 0.25;
cp_len = round(N_fft * cp_ratio);
symbol_len = N_fft + cp_len;
num_symbols = floor(t_duration * fs / symbol_len);
k_start = 2;
k_end = floor(f_end / df) + 1;
total_samples = num_symbols * symbol_len;
signal_origin_full = zeros(1, total_samples);
pos_idx = k_start : k_end;
neg_idx = N_fft - pos_idx + 2;
idx_pointer = 1;

for s = 1:num_symbols
    data_bits = randi([0, 1], length(pos_idx), 1);
    mod_data = 2 * data_bits - 1;
    
    spectrum = zeros(N_fft, 1);
    spectrum(pos_idx) = mod_data;
    spectrum(neg_idx) = conj(mod_data);
    
    ifft_sym = real(ifft(spectrum, N_fft));
    cp = ifft_sym(end - cp_len + 1 : end);
    symbol = [cp; ifft_sym];
    
    signal_origin_full(idx_pointer : idx_pointer + symbol_len - 1) = symbol.';
    idx_pointer = idx_pointer + symbol_len;
end

exact_len_full = round(t_duration * fs);
if length(signal_origin_full) > exact_len_full
    signal_origin_full = signal_origin_full(1:exact_len_full);
else
    signal_origin_full = [signal_origin_full, zeros(1, exact_len_full - length(signal_origin_full))];
end

%% --- 3. Compute Fixed Baseline Noise from 5m Tag ---
ref_dist = 5;                                      % 基准距离 5m
ref_amp = 1 / (ref_dist^4);                        % 基准衰减
avg_ook_pwr_factor = (amp_on^2 + amp_off^2) / 2;
base_sig_pwr = mean(signal_origin_full .^ 2);
p_ref = (ref_amp^2) * base_sig_pwr * avg_ook_pwr_factor;
noise_var = p_ref / (10^(snr / 10));               % 固定的全系统底噪方差
noise_std = sqrt(noise_var);

fprintf('========================================================================\n');
fprintf('  Multi-Tag System Capacity vs. Spacing & Tag Count (N = 1 to %d)\n', N_max);
fprintf('  Reference SNR at 5m: %g dB | Channel Noise Std: %e\n', snr, noise_std);
fprintf('  Target Maximum BER : <= %.4f\n', target_BER);
fprintf('========================================================================\n\n');

%% --- 4. Main Simulation Loop with Log-Space Binary Search ---
N_list = 1:N_max;
single_tag_rate_matrix = zeros(num_spacings, N_max);
total_capacity_matrix  = zeros(num_spacings, N_max);

spb_min_bound = 8;
spb_max_bound = 150000;

% 对数二分精度门限：例如 log2(1.02) 对应约 2% 的相对搜索误差
log_rel_tol = log2(1.02);
num_mc_trials = 5; % 每次判定执行 5 轮蒙特卡洛平均，抑制 BER 抽样随机性

for s_idx = 1:num_spacings
    delta_d = spacing_list(s_idx);
    fprintf('\n>>> Running Simulation for Spacing = %.1f m <<<\n', delta_d);
    
    best_spb_last = spb_min_bound;
    
    for n_idx = 1:N_max
        curr_N = N_list(n_idx);
        
        % 当前组与当前 N 下各 Tag 的位置、时延与衰减
        distances = 5 + (0 : curr_N - 1) * delta_d;
        delay_samples = round(2 * distances * fs / c);
        amplitudes = 1 ./ (distances .^ 4);
        
        % 对数域搜索边界设定
        log_low = log2(best_spb_last);
        log_high = log2(spb_max_bound);
        best_spb = -1;
        
        % --- 对数域二分搜索核心逻辑 ---
        while (log_high - log_low) > log_rel_tol
            mid_log = (log_low + log_high) / 2;
            mid_spb = round(2^mid_log);
            
            % 执行多轮平均 BER 测试
            [is_satisfied, ~] = test_multi_tag_ber_averaged(...
                signal_origin_full, mid_spb, curr_N, delay_samples, amplitudes, ...
                amp_on, amp_off, noise_std, target_BER, num_mc_trials);
            
            if is_satisfied
                best_spb = mid_spb;
                % 满足 BER 要求，说明 SPB 可以进一步减小以换取更高速率
                log_high = mid_log;
            else
                % 不满足要求，说明 SPB 太小（信噪比不够），需要加大 SPB
                log_low = mid_log;
            end
        end
        
        % 边界收敛微调：取收敛点或其邻域上限作为安全最佳解
        if best_spb < 0
            % 若在循环中未成功捕获，测试上边界
            check_spb = round(2^log_high);
            [is_sat, ~] = test_multi_tag_ber_averaged(...
                signal_origin_full, check_spb, curr_N, delay_samples, amplitudes, ...
                amp_on, amp_off, noise_std, target_BER, num_mc_trials);
            if is_sat
                best_spb = check_spb;
            end
        end
        
        if best_spb > 0
            best_spb_last = best_spb; % 动态复用作为下一个 N 的搜索下界
            r_tag = fs / best_spb;
            single_tag_rate_matrix(s_idx, n_idx) = r_tag;
            total_capacity_matrix(s_idx, n_idx)  = curr_N * r_tag;
            
            fprintf('[Spacing %.1fm | N = %2d/%2d] Far Dist = %5.1fm | SPB = %-6d | R_tag = %8.2f kbps | Cap = %8.2f kbps\n', ...
                delta_d, curr_N, N_max, distances(end), best_spb, ...
                r_tag / 1e3, total_capacity_matrix(s_idx, n_idx) / 1e3);
        else
            fprintf('[Spacing %.1fm | N = %2d/%2d] Failed to converge!\n', delta_d, curr_N, N_max);
        end
    end
end

fprintf('\n=== Simulation Completed. Generating Comparison Plots... ===\n');

%% --- 5. Visualization (IEEE Standard Comparison Format) ---
figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 750, 650]);
colors = [
    0.000, 0.447, 0.741;  % 蓝色: 0.5m
    0.850, 0.325, 0.098;  % 橙色: 0.7m
    0.466, 0.674, 0.188   % 绿色: 1.0m
];
markers = {'-o', '-s', '-^'};

% --- 子图 1: 单 Tag 速率随 N 变化 ---
subplot(2, 1, 1);
hold on;
for s_idx = 1:num_spacings
    plot(N_list, single_tag_rate_matrix(s_idx, :) / 1e3, markers{s_idx}, ...
        'LineWidth', 1.4, 'MarkerSize', 4.0, ...
        'Color', colors(s_idx, :), 'MarkerFaceColor', colors(s_idx, :), ...
        'DisplayName', sprintf('Spacing = %.1f m', spacing_list(s_idx)));
end
hold off;
xlabel('Number of Tags ($N$)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Single Tag Rate (kbps)', 'Interpreter', 'latex', 'FontSize', 11);
title('Single-Tag Modulation Rate vs. Tag Count ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9.5, 'Box', 'on');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, ...
    'LineWidth', 0.9, 'TickLabelInterpreter', 'latex', 'Box', 'on');

% --- 子图 2: 系统总容量随 N 变化 ---
subplot(2, 1, 2);
hold on;
for s_idx = 1:num_spacings
    plot(N_list, total_capacity_matrix(s_idx, :) / 1e3, markers{s_idx}, ...
        'LineWidth', 1.4, 'MarkerSize', 4.0, ...
        'Color', colors(s_idx, :), 'MarkerFaceColor', colors(s_idx, :), ...
        'DisplayName', sprintf('Spacing = %.1f m', spacing_list(s_idx)));
end
hold off;
xlabel('Number of Tags ($N$)', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Total System Capacity (kbps)', 'Interpreter', 'latex', 'FontSize', 11);
title('Total System Capacity ($N \times R_{\mathrm{tag}}$) vs. Tag Count ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 9.5, 'Box', 'on');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, ...
    'LineWidth', 0.9, 'TickLabelInterpreter', 'latex', 'Box', 'on');

%% --- Helper Function: Averaged Multi-Tag BER Check ---
function [satisfied, avg_max_ber] = test_multi_tag_ber_averaged(...
    signal_origin_full, samples_per_bit, N, delay_samples, amplitudes, ...
    amp_on, amp_off, noise_std, target_BER, num_trials)

    ber_records = zeros(1, num_trials);
    for t = 1:num_trials
        [is_sat, cur_max_ber] = test_multi_tag_ber_single(...
            signal_origin_full, samples_per_bit, N, delay_samples, amplitudes, ...
            amp_on, amp_off, noise_std, target_BER);
        
        ber_records(t) = cur_max_ber;
        % 若单次误码率已严重恶化（例如大于 5 倍门限），可提前熔断以加速仿真
        if cur_max_ber > 5 * target_BER
            satisfied = false;
            avg_max_ber = cur_max_ber;
            return;
        end
    end
    
    avg_max_ber = mean(ber_records);
    satisfied = (avg_max_ber <= target_BER);
end

%% --- Helper Function: Exact SIC Demodulation Single Trial ---
function [satisfied, max_ber] = test_multi_tag_ber_single(...
    signal_origin_full, samples_per_bit, N, delay_samples, amplitudes, ...
    amp_on, amp_off, noise_std, target_BER)

    max_delay = max(delay_samples);
    num_bits = floor(length(signal_origin_full) / samples_per_bit);
    
    if num_bits < 30
        satisfied = false;
        max_ber = 1.0;
        return;
    end
    
    exact_len = num_bits * samples_per_bit;
    signal_origin = signal_origin_full(1:exact_len);
    
    % 生成随机 bit 与对应幅度序列
    bits = randi([0, 1], N, num_bits);
    amp_seq = bits * (amp_on - amp_off) + amp_off;
    
    total_len = exact_len - max_delay;
    rx_mix = zeros(1, total_len);
    
    % 多 Tag 回波在时域线性叠加
    for i = 1:N
        tau = delay_samples(i);
        ook_sym = repelem(amp_seq(i, :), 1, samples_per_bit);
        delayed_ref = signal_origin(tau + 1 : total_len + tau);
        rx_mix = rx_mix + amplitudes(i) * delayed_ref .* ook_sym(1:total_len);
    end
    
    % 叠加恒定绝对基准白噪声
    noise = noise_std * randn(1, total_len);
    rx_mix = rx_mix + noise;
    
    % --- 精准递归 SIC 解调 (按信号由强到弱排序) ---
    [~, sic_order] = sort(amplitudes, 'descend');
    
    demod_bits = zeros(N, num_bits);
    current_rx = rx_mix;
    num_chunks = floor(total_len / samples_per_bit);
    
    for rank_idx = 1:N
        tag_id = sic_order(rank_idx);
        tau = delay_samples(tag_id);
        amp_k = amplitudes(tag_id);
        
        threshold = amp_k * (amp_on + amp_off) / 2;
        reconstructed_ook = zeros(1, num_bits);
        
        for b = 1:num_chunks
            idx_start = (b - 1) * samples_per_bit + 1;
            idx_end = b * samples_per_bit;
            
            rx_seg = current_rx(idx_start : idx_end);
            ref_seg = signal_origin(tau + idx_start : tau + idx_end);
            
            ref_energy = sum(ref_seg .^ 2);
            if ref_energy > 0
                dot_val = sum(rx_seg .* ref_seg) / ref_energy;
            else
                dot_val = 0;
            end
            
            % OOK 硬判决
            if dot_val >= threshold
                demod_bits(tag_id, b) = 1;
                reconstructed_ook(b) = amp_on;
            else
                demod_bits(tag_id, b) = 0;
                reconstructed_ook(b) = amp_off;
            end
        end
        
        % SIC: 重构并消去当前 Tag 回波
        reconstructed_full_ook = repelem(reconstructed_ook, 1, samples_per_bit);
        reconstructed_sub_sig = amp_k * signal_origin(tau + 1 : total_len + tau) .* reconstructed_full_ook(1:total_len);
        current_rx = current_rx - reconstructed_sub_sig;
    end
    
    % 统计全体 Tag 误码率
    ber_per_tag = zeros(1, N);
    for i = 1:N
        errs = sum(demod_bits(i, 1:num_chunks) ~= bits(i, 1:num_chunks));
        ber_per_tag(i) = errs / num_chunks;
    end
    
    max_ber = max(ber_per_tag);
    satisfied = (max_ber <= target_BER);
end