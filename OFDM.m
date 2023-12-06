% ͨ��ϵͳ������SOC���� 2023���＾ѧ��
clear;clc;close all;

% ��ȡ��Ч�����ŵ��弤��Ӧh
load('h.mat');

% SNR��λΪdB�������е���
SNR_list = 0 : 2 : 20;

% ÿ��SNR���������OFDM���������������е���
block_num = 10000;

%��Ƶ���
pilot_interval = 5;

% FFT����
N_fft = 1024;     

% ���ƽ�����4(QPSK) �� 16(16QAM)
Q = 16;

% ÿ�����ų��صı�����
B = log2(Q);

% ����ؼ���
BER_count = zeros(size(SNR_list));

% ������Ϣ�����ز�����,�趨Ϊ�������������ز���
N_sc = 716;

% ÿ��OFDM����������ı�����
N_bit = N_sc * B;

% CP����
length_CP = 73;

% equalization method, 0: ZF, 1:MMSE
equal_method = 1;

for snr_count = 1:length(SNR_list)
    
    snr = SNR_list(snr_count);
    
    % ���������������Ƶ�������ƽ�����ʹ�һ��Ϊ1��SNR�Ķ���Ϊ��N_sc / (N_fft ^ 2 * ��������)
    sigma2 = N_sc / (10 ^ (snr / 10) * N_fft ^ 2);
    
    for count = 1 : block_num
        snr_count, count
        
        %������Ϣ����
        msg = round(rand(1, N_bit));
        
        %% OFDM���Ʋ���
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬��ɴӱ��ص�QPSK��16QAM���������ӳ�䣬���뵼Ƶ���ţ�����ĵ�Ƶ���Ź���Ϊ1��������OFDM���Ʋ����CP
        % ���룺��Ϣ��������msg
        % �����OFDM���ƺ����CP��ʱ���ź�����x_ofdm������������ʾ
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

        %% �ŵ����䲿��
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % �������Ĵ��룬�������x_ofdm�����ྶ�ŵ�h������նˣ������AWGN��ȥ��CP��ע��ʵ�鲿�������ʸ�Ϊ���������ʵ�һ��
        % �����Ѿ������CP����������һ��OFDM���ŶԱ�OFDM���ŵ�Ӱ��
        % ���룺���CP���ʱ��OFDM�������x_ofdm
        % ���������ΪN_fft�Ľ���OFDM����r_ofdm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        x_h = conv(x_ofdm, h);
        x_h_awgn = x_h + randn(1, length(x_h)) * sqrt(0.5 * sigma2) + randn(1, length(x_h)) * sqrt(0.5 * sigma2) * 1i;
        r_ofdm = x_h_awgn(length_CP+1:length_CP+N_fft);
        
        %% OFDM�������       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ��������Լ��Ĵ��룬���õ�Ƶ�����ŵ�������ʹ�ù��Ƶõ����ŵ��Խ��յ���OFDM���Ž��������MMSEƵ����⣬�о������㲢�ָ��������
        % ���룺����OFDM����r_ofdm
        % ����������ı�������msg_r
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

        %% �������ͳ��
        BER_count(snr_count) = BER_count(snr_count) + sum(abs(msg_r - msg));
        
    end
    
end

% �����ʼ���
BER = BER_count / (block_num * N_bit);

%% BER���߻���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������Ĵ��룬���ð�������꣬ʹ��semilogy��������BER-SNR����
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
