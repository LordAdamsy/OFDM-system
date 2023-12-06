% 通信系统仿真与SOC集成 2023年秋季学期
clear;clc;close all;

% 读取等效基带信道冲激响应h
load('h.mat');

% SNR单位为dB，可自行调整
SNR_list = 0 : 2 : 20;

% 每个SNR下所传输的OFDM符号数量，可自行调整
block_num = 10000;

%导频间隔
pilot_interval = 5;

% FFT点数
N_fft = 1024;     

% 调制阶数：4(QPSK) 或 16(16QAM)
Q = 16;

% 每个符号承载的比特数
B = log2(Q);

% 误比特计数
BER_count = zeros(size(SNR_list));

% 传输信息的子载波数量,设定为计算出的最大子载波数
N_sc = 716;

% 每个OFDM符号所传输的比特数
N_bit = N_sc * B;

% CP点数
length_CP = 73;

% equalization method, 0: ZF, 1:MMSE
equal_method = 1;

for snr_count = 1:length(SNR_list)
    
    snr = SNR_list(snr_count);
    
    % 计算噪声方差，其中频域发射符号平均功率归一化为1，SNR的定义为：N_sc / (N_fft ^ 2 * 噪声方差)
    sigma2 = N_sc / (10 ^ (snr / 10) * N_fft ^ 2);
    
    for count = 1 : block_num
        snr_count, count
        
        %生成信息比特
        msg = round(rand(1, N_bit));
        
        %% OFDM调制部分
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你的代码，完成从比特到QPSK（16QAM）星座点的映射，插入导频符号（插入的导频符号功率为1），进行OFDM调制并添加CP
        % 输入：信息比特序列msg
        % 输出：OFDM调制后并添加CP的时域信号序列x_ofdm，用行向量表示
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % preparation for modulation
        msg_mid = reshape(msg, B, []);

        % qpsk modulation or 16qam modulation
        if Q == 4
            msg_sym_int = bit2int(msg_mid, B);
            msg_mod = pskmod(msg_sym_int, Q);
        elseif Q == 16
            msg_sym_int = bit2int(msg_mid, B);
            msg_mod = 1 /sqrt(10) * qammod(msg_sym_int, Q);            
        end
        
        % S to D
        row = N_fft;
        col = ceil(length(msg_mod)/row);

        % data matrix
        data = zeros(row, col);

        % compute index of pilot sequence and data
        idx = 1:N_fft;
        pilot_idx = 1:pilot_interval:N_fft;
        data_idx = idx(~ismember(idx, pilot_idx));

        % input data and pilot matrix, P_pilot = 1
        pilot_matrix = (1+1i) / sqrt(2) * ones(length(pilot_idx), col);
        data(pilot_idx, :) = pilot_matrix;
        msg_mod_vec = [msg_mod, zeros(1, (N_fft-length(pilot_idx))*col-length(msg_mod))];
        msg_mod_mat = reshape(msg_mod_vec, [], col);
        data(data_idx) = msg_mod_mat;

        % ofdm mod by ifft
        data_ifft = ifft(data);
        
        % D to S
        data_ifft_seq = reshape(data_ifft, 1, []);

        %cyclic prefix
        x_ofdm = [data_ifft_seq(length(data_ifft_seq)-length_CP+1:end), data_ifft_seq];

        %% 信道传输部分
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你的代码，发射符号x_ofdm经过多径信道h到达接收端，并添加AWGN、去掉CP。注意实虚部噪声功率各为总噪声功率的一半
        % 由于已经添加了CP，不考虑上一个OFDM符号对本OFDM符号的影响
        % 输入：添加CP后的时域OFDM发射符号x_ofdm
        % 输出：长度为N_fft的接收OFDM符号r_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        x_h = conv(x_ofdm, h);
        x_h_awgn = x_h + randn(1, length(x_h)) * sqrt(0.5 * sigma2) + randn(1, length(x_h)) * sqrt(0.5 * sigma2) * 1i;
        r_ofdm = x_h_awgn(length_CP+1:length_CP+N_fft);
        
        %% OFDM解调部分       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 请添加你自己的代码，利用导频估计信道，进而使用估计得到的信道对接收到的OFDM符号进行迫零或MMSE频域均衡，判决星座点并恢复传输比特
        % 输入：接收OFDM符号r_ofdm
        % 输出：解调后的比特序列msg_r
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % S to D
        r_ofdm_d = reshape(r_ofdm, N_fft, []);

        % fft
        r_fft = fft(r_ofdm_d);

        % channel estimate
        pilot_matrix_r = r_fft(pilot_idx, :);
        h_pilot = pilot_matrix_r ./ pilot_matrix;
        h_data = interp1(pilot_idx, h_pilot, data_idx).';

        % equalization
        if equal_method == 0 % ZF equalization
            x_equal = r_fft(data_idx, :) ./ h_data;
        elseif equal_method == 1 % MMSE equalization
            x_equal = r_fft(data_idx, :) .* conj(h_data) ./ (abs(h_data) .^ 2 + 1/(10 ^ (snr / 10)));
        end

        % D to S
        data_r_zeros = reshape(x_equal, 1, []);

        % delete zeros
        data_r = data_r_zeros(1:N_sc);

        % demodulation
        if Q == 4
            msg_r_int = pskdemod(data_r, Q).';
        elseif Q == 16
            msg_r_int = qamdemod(sqrt(10) * data_r, Q).';
        end

        msg_r_bit = int2bit(msg_r_int, B);
        msg_r = reshape(msg_r_bit, 1, []);

        %% 误比特数统计
        BER_count(snr_count) = BER_count(snr_count) + sum(abs(msg_r - msg));
        
    end
    
end

% 误码率计算
BER = BER_count / (block_num * N_bit);

%% BER曲线绘制
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 请添加你的代码，采用半对数坐标，使用semilogy函数绘制BER-SNR曲线
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Q == 4 && equal_method == 0
    figure
    semilogy(SNR_list, BER, 'LineWidth', 1.5);
    xlabel('SNR/dB');
    ylabel('BER'); 
    title('BER of OFDM system with QPSK modulation and ZF equlization');
    grid on;
elseif Q == 4 && equal_method == 1
    figure
    semilogy(SNR_list, BER, 'LineWidth', 1.5);
    xlabel('SNR/dB');
    ylabel('BER');
    title('BER of OFDM system with QPSK modulation and MMSE equlization');
    grid on;
elseif Q == 16 && equal_method == 0
    figure
    semilogy(SNR_list, BER, 'LineWidth', 1.5);
    xlabel('SNR/dB');
    ylabel('BER');
    title('BER of OFDM system with 16QAM modulation and ZF equlization');
    grid on;
else
    figure
    semilogy(SNR_list, BER, 'LineWidth', 1.5);
    xlabel('SNR/dB');
    ylabel('BER');
    title('BER of OFDM system with 16QAM modulation and MMSE equlization');
    grid on;
end
