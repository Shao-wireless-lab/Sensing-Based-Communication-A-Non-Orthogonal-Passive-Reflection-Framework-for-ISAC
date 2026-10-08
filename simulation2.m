% Comprehensive Simulation: BER and Localization Error vs. Bit Rate (With Known Location Baseline)
clc; close all; clear;

%% --- 1. Parameter Settings ---
fs = 2.5e9;                   % System sampling rate (2.5 GHz)
t_duration = 0.01;            % Signal duration (s)
f_start = 10000;              % Start frequency (10 kHz)
f_end = 250e6;                % End frequency (250 MHz)
bandwidth = f_end - f_start;  % Signal bandwidth

%% --- 2. OFDM & Subcarrier Design ---
N_fft = 1024;                     % FFT size
df = fs / N_fft;                  % Subcarrier spacing
cp_ratio = 0.25;                  % Cyclic prefix (CP) ratio
cp_len = round(N_fft * cp_ratio); % CP length
symbol_len = N_fft + cp_len;      % Total OFDM symbol length
num_symbols = floor(t_duration * fs / symbol_len); 
k_start = 2;                          
k_end = floor(f_end / df) + 1;        
active_subcarriers = k_end - k_start + 1; 

%% --- 3. Transmitted Signal Generation (Base OFDM Waveform) ---
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
    
    ifft_sym = ifft(spectrum, N_fft);
    ifft_sym = real(ifft_sym); 
    
    cp = ifft_sym(end - cp_len + 1 : end);
    symbol = [cp; ifft_sym];
    
    signal_origin_full(idx_pointer : idx_pointer + symbol_len - 1) = symbol.';
    idx_pointer = idx_pointer + symbol_len;
end
exact_len_full = round(t_duration * fs);
if length(signal_origin_full) > exact_len_full
    signal_origin_full = signal_origin_full(1:exact_len_full);
elseif length(signal_origin_full) < exact_len_full
    signal_origin_full = [signal_origin_full, zeros(1, exact_len_full - length(signal_origin_full))];
end

%% --- 4. Bit Rate and Samples-Per-Bit Grid ---
max_samples = 200;
min_samples = 8;
min_bit_rate = fs / max_samples; 
max_bit_rate = fs / min_samples;  
target_bit_rates = linspace(min_bit_rate, max_bit_rate, 60); 
samples_list = round(fs ./ target_bit_rates); 
samples_list = unique(samples_list, 'stable'); 
num_points = length(samples_list);

% Preallocate metric arrays
bit_rate_list              = zeros(1, num_points);
p95_loc_error_list         = zeros(1, num_points);
median_loc_error_list      = zeros(1, num_points);
ber_list                   = zeros(1, num_points);
ber_known_loc_list         = zeros(1, num_points); % 新增：位置已知时的 BER

% Physical channel parameters
c = 3e8;                                      % Speed of light (m/s)
distance = 5;                                 % Distance (m)
delay_samples = round(2 * distance * fs / c); % Round-trip sample delay
amplitude = 1 / (distance^4); 
snr = 15;                                     % Channel SNR (dB)
cm_per_sample = (c / (2 * fs)) * 100;         % One-way distance per sample (cm)
fprintf('=== Starting Simulation (%d test points) ===\n', num_points);

%% --- 5. Main Simulation Loop ---
for iter = 1:num_points
    samples_per_bit = samples_list(iter);
    bit_rate = fs / samples_per_bit;
    bit_rate_list(iter) = bit_rate;
    
    fprintf('[%2d/%2d] Processing samples_per_bit = %-3d | Bit Rate = %-6.2f Mbps\n', ...
            iter, num_points, samples_per_bit, bit_rate / 1e6);
            
    amp_on = 1;
    amp_off = 0.1;
    num_tags = 1;
    
    num_bits = floor(exact_len_full / samples_per_bit); 
    exact_len = num_bits * samples_per_bit;
    
    signal_origin = signal_origin_full(1:exact_len); 
    bits = randi([0, 1], num_tags, num_bits); 
    amp_sequence = bits * (amp_on - amp_off) + amp_off;
    
    % OOK symbol repetition
    ook_signal = repelem(amp_sequence, 1, samples_per_bit);
    signal = signal_origin;
    
    % Propagation & Backscatter
    total_length = length(signal) - 1000;
    rx_signal = zeros(1, total_length);
    rx_signal(1, :) = amplitude * signal(1, delay_samples + 1 : total_length + delay_samples);
    
    len_rx = length(rx_signal(1, :)); 
    rx_signal(1, :) = rx_signal(1, :) .* ook_signal(1:len_rx);
    rx_signal(1, :) = awgn(rx_signal(1, :), snr, 'measured');
    
    % Demodulation & Sliding Cross-Correlation
    mix_signal = rx_signal; 
    limit = 200;
    len_mix = samples_per_bit;
    num_chunks = floor((total_length - limit) / len_mix);
    
    estimated_delays      = zeros(1, num_chunks);
    demodulated_ook       = zeros(1, num_chunks);
    demodulated_ook_known = zeros(1, num_chunks); % 新增：已知位置下的判决序列
    decision_threshold    = amplitude * (amp_on + amp_off) / 2;
    
    for c_idx = 1:num_chunks
        start_idx = (c_idx - 1) * len_mix + 1;
        mix_signal_segment = mix_signal(start_idx : start_idx + len_mix - 1);
        
        % --- 动态估计位置路径 ---
        segment_sig = signal_origin(start_idx : start_idx + len_mix + limit - 2);
        corr_results = conv(segment_sig, fliplr(mix_signal_segment), 'valid');
        
        [raw_max_val, max_idx] = max(corr_results);
        
        best_ref_start = start_idx + max_idx - 1;
        best_signal_segment = signal_origin(best_ref_start : best_ref_start + len_mix - 1);
        ref_energy = sum(best_signal_segment .^ 2);
        
        if ref_energy > 0
            normalized_max = raw_max_val / ref_energy;
        else
            normalized_max = 0;
        end
        
        estimated_delays(c_idx) = max_idx - 1; 
        demodulated_ook(c_idx) = (normalized_max >= decision_threshold);
        
        % --- 物体位置已知路径（无估计抖动，直接采用准确延迟采样点）---
        known_ref_start = start_idx + delay_samples;
        known_ref_segment = signal_origin(known_ref_start : known_ref_start + len_mix - 1);
        known_ref_energy = sum(known_ref_segment .^ 2);
        
        % 在精确延迟位置计算点积相关值
        corr_known = sum(mix_signal_segment .* known_ref_segment);
        if known_ref_energy > 0
            normalized_known = corr_known / known_ref_energy;
        else
            normalized_known = 0;
        end
        demodulated_ook_known(c_idx) = (normalized_known >= decision_threshold);
    end
    
    % --- Metric Calculation ---
    true_bits = bits(1, 1:num_chunks);
    
    % 1. Bit Error Rate (BER)
    ber_list(iter) = sum(demodulated_ook ~= true_bits) / num_chunks;
    ber_known_loc_list(iter) = sum(demodulated_ook_known ~= true_bits) / num_chunks;
    
    % 2. Localization Error (95th Percentile & Median)
    delay_error_samples = abs(estimated_delays - delay_samples);
    loc_errors_cm = delay_error_samples * cm_per_sample;
    
    sorted_errors = sort(loc_errors_cm);
    idx_95 = max(1, round(0.95 * length(sorted_errors)));
    idx_50 = max(1, round(0.50 * length(sorted_errors)));
    
    p95_loc_error_list(iter) = sorted_errors(idx_95);
    median_loc_error_list(iter) = sorted_errors(idx_50);
end
fprintf('\n=== Simulation Completed. ===\n');

%% --- 5.5 理论计算与静态参数应用 (直接代入确定的 P_peak_th) ---
fprintf('=== Calculating Theoretical Model (Direct Implementation) ===\n');

% 基础理论参数 (与仿真保持严格一致)
B_th = 250e6 - 10000;       
beta_th = 0.1;              
A_n_th = 0.04 * 1/(5^4);    
d_th = 5;                   
a_th = 1 / (d_th^4);           
Interf_th = A_n_th^2; 

% 按比特率升序排列仿真结果 (用于对齐画图用的 X 轴范围)
[bit_rate_sorted, sort_idx] = sort(bit_rate_list);
ber_known_loc_sorted = ber_known_loc_list(sort_idx);
ber_sorted = ber_list(sort_idx);

% 直接指定已知参数，不再进行优化搜索
P_peak_known = 0.027335;
P_peak_dyn   = 0.404890;

% 定义匿名函数1：根据 P_peak_th 和速率 Rs 计算理论 BER (保留原始推导公式)
calc_theo_ber = @(P_peak, Rs) 0.5 * erfc( sqrt(B_th ./ Rs) .* (a_th * (1 - beta_th)) ./ ...
    ( sqrt(2 * P_peak * a_th^2 + Interf_th) + sqrt(2 * (beta_th^2) * P_peak * a_th^2 + Interf_th) ) / sqrt(2) );

% 利用确定的 P_peak_th 生成平滑的理论曲线以供画图
Rs_range_th = linspace(min(bit_rate_sorted), max(bit_rate_sorted), 500); 
ber_th_known_vals = calc_theo_ber(P_peak_known, Rs_range_th);
ber_th_dyn_vals   = calc_theo_ber(P_peak_dyn, Rs_range_th);

%% --- 6. Plotting (IEEE Standard Format) ---
p95_loc_error_sorted    = p95_loc_error_list(sort_idx);
median_loc_error_sorted = median_loc_error_list(sort_idx);
bit_rate_mbps           = bit_rate_sorted / 1e6;

figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 700, 600]);

% --- Subplot 1: Localization Error vs. Bit Rate ---
subplot(2, 1, 1);
plot(bit_rate_mbps, p95_loc_error_sorted, '-o', ...
    'LineWidth', 1.4, 'MarkerSize', 4.0, ...
    'Color', [0.0, 0.447, 0.741], 'MarkerFaceColor', [0.0, 0.447, 0.741]);
hold on;
plot(bit_rate_mbps, median_loc_error_sorted, '--s', ...
    'LineWidth', 1.2, 'MarkerSize', 4.0, ...
    'Color', [0.85, 0.325, 0.098], 'MarkerFaceColor', [0.85, 0.325, 0.098]);
hold off;
xlabel('Bit Rate (Mbps)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Localization Error (cm)', 'Interpreter', 'latex', 'FontSize', 18);
legend({'Dynamic 95th Percentile Error', 'Dynamic Median Error (50th)'}, ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'northwest', 'Box', 'on');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, ...
    'LineWidth', 0.9, 'TickLabelInterpreter', 'latex', 'Box', 'on');

% --- Subplot 2: Bit Error Rate (BER) vs. Bit Rate ---
subplot(2, 1, 2);
ber_floor = 1e-6; % 对数坐标下零误码的展示下限
ber_display = ber_sorted;
ber_display(ber_display == 0) = ber_floor;
ber_known_display = ber_known_loc_sorted;
ber_known_display(ber_known_display == 0) = ber_floor;

% 1. 绘制仿真：动态估计 BER
semilogy(bit_rate_mbps, ber_display, '-^', ...
    'LineWidth', 1.4, 'MarkerSize', 4.0, ...
    'Color', [0.635, 0.078, 0.184], 'MarkerFaceColor', [0.635, 0.078, 0.184]);
hold on;
% 2. 绘制仿真：理想同步(Known Location) BER
semilogy(bit_rate_mbps, ber_known_display, '--d', ...
    'LineWidth', 1.4, 'MarkerSize', 4.0, ...
    'Color', [0.466, 0.674, 0.188], 'MarkerFaceColor', [0.466, 0.674, 0.188]);

% 3. 绘制理论拟合曲线：动态估计
semilogy(Rs_range_th / 1e6, ber_th_dyn_vals, '-.k', 'LineWidth', 2.0);
% 4. 绘制理论拟合曲线：已知位置
semilogy(Rs_range_th / 1e6, ber_th_known_vals, ':k', 'LineWidth', 2.0);

hold off;
xlabel('Bit Rate (Mbps)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', 18);

% 简化的图例，不再标出 P_peak 的具体数值
legend({'Dynamic Estimation (Sim)', 'Known Location (Sim)', 'Theory Dynamic', 'Theory Known'}, ...
    'Interpreter', 'latex', 'FontSize', 15, 'Location', 'southwest', 'Box', 'on');
ylim([ber_floor, 1]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, ...
    'LineWidth', 0.9, 'TickLabelInterpreter', 'latex', 'Box', 'on');