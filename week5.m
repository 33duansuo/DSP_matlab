clear
clc
close all
%% 程序版本信息
% 作者： 段仁俊
% 时间： 20241130
%% 设定数据文件保存地址
FolderPath='E:\课程\DSP\';
% 数据文件名
Save_Path = strcat(FolderPath,'Result2\');
% 处理结果保存地址
%% ===================【程序说明】===================
% 无
%% 第一题
% 参数设置
close all
N = 45; % 窗函数的长度

% 计算不同窗函数
rect_window = ones(1, N); % 矩形窗
hamming_window = hamming(N)'; % 汉明窗
blackman_window = blackman(N)'; % 布莱克曼窗

% 使用 freqz 计算每个窗函数的频率响应
[rect_h, rect_f] = freqz(rect_window, 1); % 矩形窗频率响应
[hamming_h, hamming_f] = freqz(hamming_window, 1); % 汉明窗频率响应
[blackman_h, blackman_f] = freqz(blackman_window, 1); % 布莱克曼窗频率响应
% 计算幅度谱
rect_mag = 20*log10(abs(rect_h)/ max(abs(rect_mag))); % 矩形窗幅度谱
hamming_mag = 20*log10(abs(hamming_h)/ max(abs(rect_mag))); % 汉明窗幅度谱
blackman_mag =20*log10(abs(blackman_h)/ max(abs(rect_mag))); % 布莱克曼窗幅度谱

% 绘制频谱
figure;

subplot(3, 1, 1);
plot(rect_f/pi, rect_mag);
title('矩形窗的归一化幅度谱');
xlabel('归一化频率\pi');
ylabel('幅度/dB');
grid on;

subplot(3, 1, 2);
plot(hamming_f/pi, hamming_mag);
title('汉明窗的归一化幅度谱');
xlabel('归一化频率\pi');
ylabel('幅度/dB');
grid on;

subplot(3, 1, 3);
plot(2*blackman_f/pi, blackman_mag);
title('布莱克曼窗的归一化幅度谱');
xlabel('归一化频率\pi');
ylabel('幅度/dB');
grid on;
%% 保存图片
saveas(1, [Save_Path '第一问'],'svg'); 
close all;
%% 第二题
% 参数设置
close all
% 参数设置
N1 = 15; % 滤波器长度 N = 15
N2 = 45; % 滤波器长度 N = 45
omega1 = 0.3 * pi; % 通带下边界
omega2 = 0.5 * pi; % 通带上边界
fs = 1; % 归一化频率
f = (0:1023)/1024 * fs; % 频率轴

% 使用 fir1 和汉宁窗设计带通滤波器 (N=15)
h_15 = fir1(N1-1, [omega1/pi omega2/pi], 'bandpass', hanning(N1)); % 带通滤波器
[H_15, F_15] = freqz(h_15, 1, 1024, fs); % 计算滤波器的频率响应

% 使用 fir1 和汉宁窗设计带通滤波器 (N=45)
h_45 = fir1(N2-1, [omega1/pi omega2/pi], 'bandpass', hanning(N2)); % 带通滤波器
[H_45, F_45] = freqz(h_45, 1, 1024, fs); % 计算滤波器的频率响应

% 绘制频率响应
figure;

subplot(2, 2, 1);
plot(F_15*2, 20*log10(abs(H_15)));
title('N = 15，带通滤波器的幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;

subplot(2, 2, 2);
plot(F_15*2, angle(H_15));
title('N = 15，带通滤波器的相位响应');
xlabel('归一化频率\pi');
ylabel('相位 (radians)');
grid on;

subplot(2, 2, 3);
plot(F_45*2, 20*log10(abs(H_45)));
title('N = 45，带通滤波器的幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;

subplot(2, 2, 4);
plot(F_45*2, angle(H_45));
title('N = 45，带通滤波器的相位响应');
xlabel('归一化频率\pi');
ylabel('相位 (radians)');
grid on;

% 计算和显示3dB和20dB带宽
[bw_15_3dB, bw_15_20dB] = calculate_bandwidth(H_15, F_15);
[bw_45_3dB, bw_45_20dB] = calculate_bandwidth(H_45, F_45);

fprintf('N=15时的3dB带宽: %.4f\n', bw_15_3dB);
fprintf('N=15时的20dB带宽: %.4f\n', bw_15_20dB);
fprintf('N=45时的3dB带宽: %.4f\n', bw_45_3dB);
fprintf('N=45时的20dB带宽: %.4f\n', bw_45_20dB);
%% 保存图片
saveas(1, [Save_Path '第二问'],'svg'); 
close all
%% 第三题
% 参数设置
close all
N1 = 15; % 滤波器长度 N = 15
N2 = 45; % 滤波器长度 N = 45
omega1 = 0.3 * pi; % 通带下边界
omega2 = 0.5 * pi; % 通带上边界

% 使用矩形窗设计带通滤波器 (N=15)
h_rect_15 = fir1(N1-1, [omega1/pi omega2/pi], 'bandpass', rectwin(N1)); % 矩形窗
[H_rect_15, F_rect_15] = freqz(h_rect_15, 1); % 计算频率响应

% 使用布莱克曼窗设计带通滤波器 (N=15)
h_blackman_15 = fir1(N1-1, [omega1/pi omega2/pi], 'bandpass', blackman(N1)); % 布莱克曼窗
[H_blackman_15, F_blackman_15] = freqz(h_blackman_15, 1); % 计算频率响应

% 使用矩形窗设计带通滤波器 (N=45)
h_rect_45 = fir1(N2-1, [omega1/pi omega2/pi], 'bandpass', rectwin(N2)); % 矩形窗
[H_rect_45, F_rect_45] = freqz(h_rect_45, 1); % 计算频率响应

% 使用布莱克曼窗设计带通滤波器 (N=45)
h_blackman_45 = fir1(N2-1, [omega1/pi omega2/pi], 'bandpass', blackman(N2)); % 布莱克曼窗
[H_blackman_45, F_blackman_45] = freqz(h_blackman_45, 1); % 计算频率响应

% 绘制频率响应
figure;

% N=15 矩形窗
subplot(2, 2, 1);
plot(F_rect_15*2, 20*log10(abs(H_rect_15)));
title('N=15 矩形窗幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;

% N=15 布莱克曼窗
subplot(2, 2, 2);
plot(F_blackman_15*2, 20*log10(abs(H_blackman_15)));
title('N=15 布莱克曼窗幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;

% N=45 矩形窗
subplot(2, 2, 3);
plot(F_rect_45*2, 20*log10(abs(H_rect_45)));
title('N=45 矩形窗幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;

% N=45 布莱克曼窗
subplot(2, 2, 4);
plot(F_blackman_45*2, 20*log10(abs(H_blackman_45)));
title('N=45 布莱克曼窗幅度响应');
xlabel('归一化频率\pi');
ylabel('幅度 (dB)');
grid on;
%% 保存图片
saveas(1, [Save_Path '第三问'],'svg'); 
close all
%% 第四题
close all
% 参数设置
N = 40; % 滤波器长度
beta_values = [4, 6, 10]; % 凯塞窗的不同beta值
% 带通滤波器的频率边界（归一化频率，0到1对应0到pi）
bands = [0.2 0.4 0.6 0.8]; % 归一化频率：0.2π到0.4π和0.6π到0.8π
% 创建一个图形窗口用于绘图
figure;

for i = 1:length(beta_values)
    beta = beta_values(i);
    
    % 生成凯塞窗
    window = kaiser(N, beta);
    
    % 使用fir1设计带通滤波器，'bandpass'指定带通滤波器
    b = fir1(N-1, bands, 'bandpass', window); 
    
    % 计算频率响应
    [H, W] = freqz(b, 1, 1024, fs);
    H_mag = 20*log10(abs(H));
    H_phase = angle(H);
    
    % 绘制幅度响应
    subplot(length(beta_values), 2, 2*i-1);
    plot(2*W, H_mag, 'LineWidth', 1.5);
    title(['N=40, Kaiser Window (\beta = ', num2str(beta), ') - 幅度响应']);
    xlabel('归一化频率 (\times\pi rad/sample)');
    ylabel('幅度 (dB)');
    grid on;
    axis([0 1 -100 5]); % 设置幅度轴范围
    
    % 绘制相位响应
    subplot(length(beta_values), 2, 2*i);
    plot(W*2, H_phase, 'LineWidth', 1.5);
    title(['N=40, Kaiser Window (\beta = ', num2str(beta), ') - 相位响应']);
    xlabel('归一化频率 (\times\pi rad/sample)');
    ylabel('相位 (弧度)');
    grid on;
    axis([0 1 -pi pi]); % 设置相位轴范围
end
%% 第五题
clear all
clc
close all
% 参数设置
N = 40; % 滤波器长度

% 设置幅度响应
H = zeros(1, N); % 初始化频率响应数组

% 设定通带范围的索引
k_1 = 4;   % 对应0.2π
k_2 = 8;   % 对应0.4π
k_3 = 12;  % 对应0.6π
k_4 = 16;  % 对应0.8π

% 设置幅度响应
H(k_1:k_2) = 1;  % 通带1 0.2π 到 0.4π
H(k_3:k_4) = 1;  % 通带2 0.6π 到 0.8π

% 设置过渡带范围中的幅度响应
H(k_2 + 1) = 0.5; % 过渡带1, 0.4π 到 0.5π 边界点，幅度为0.5
H(k_3 - 1) = 0.5; % 过渡带2, 0.5π 到 0.6π 边界点，幅度为0.5

% 设置阻带范围中的幅度响应
% 其余部分为0，已通过初始化为0

% 设置相位响应 θk，确保线性相位
theta = -pi * (0:N-1) * (1 - 1/N);  % 线性相位响应

% 构造频率响应 H(k)
H_k = H .* exp(1i * theta);  % 幅度与相位结合

% 计算冲激响应（通过逆离散傅里叶变换）
h = real(ifft(H_k));  % 取实部

% 构造偶对称冲激响应
% 偶对称冲激响应：h[n] = h[N-n]
h(21:end) = flip(h(1:20));  % 直接将前半部分对称到后半部分

% 绘制冲激响应
figure;
subplot(2,1,1);
stem(0:N-1, h, 'filled');
title('第二类线性相位带通滤波器的冲激响应');
xlabel('n');
ylabel('h[n]');
grid on;

% 计算频率响应并绘制幅频和相频特性
[H_freq, W] = freqz(h, 1);  % 计算频率响应

% 绘制幅频特性
subplot(2,1,2);
plot(W/pi, 20*log10(abs(H_freq)), 'LineWidth', 1.5);
title('第二类线性相位带通滤波器的幅度响应');
xlabel('归一化频率 (\times\pi rad/sample)');
ylabel('幅度 (dB)');
axis([0 1 -100 5]); % 设置幅度轴范围
grid on;

% 绘制相位特性
figure;
plot(W/pi, (angle(H_freq)), 'LineWidth', 1.5);
title('第二类线性相位带通滤波器的相位响应');
xlabel('归一化频率 (\times\pi rad/sample)');
ylabel('相位 (弧度)');
xlim([0, 1]); % 设置x轴范围从0到π
grid on;

%% 保存图片
saveas(1, [Save_Path '第五问1'],'svg'); 
saveas(2, [Save_Path '第五问2'],'svg'); 
close all
%%
clear all
close all
clc
% 设计参数
N = 40;  % 滤波器长度
f11=0.15;
f1 = 0.2; % 第一个通带起始频率（单位为π）
f2 = 0.4; % 第一个通带结束频率（单位为π）
f22=0.45;
f33=0.55;
f3 = 0.6; % 第二个通带起始频率（单位为π）
f4 = 0.8; % 第二个通带结束频率（单位为π）
f44=0.85;

% 频率采样点
frequencies = [0 f11 f1 f2 f22 f33 f3 f4 f44 1];  % 频率区间 (单位为π)
amplitude = [0 0 1 1 0  0 1 1 0 0];        % 幅度响应

% 使用firpm设计带通滤波器
b = firpm(N-1, frequencies, amplitude);  % firpm设计滤波器

% 绘制幅频特性
[H, f] = freqz(b, 1);  % 计算频率响应
f = f/pi;  % 转换为单位π
figure;
subplot(2,1,1);
plot(f, 20*log10(abs(H)));  % 幅度谱
title('带通滤波器幅度响应');
xlabel('归一化频率 (\times\pi)');
ylabel('幅度');

% 绘制相频特性
subplot(2,1,2);
plot(f, angle(H));  % 相位响应
title('带通滤波器相位响应');
xlabel('归一化频率 (\times\pi)');
ylabel('相位 (radians)');
%% 保存图片
FolderPath='E:\课程\DSP\';
% 数据文件名
Save_Path = strcat(FolderPath,'Result2\');
% 处理结果保存地址
saveas(1, [Save_Path '第六问'],'svg'); 
close all
%% 第七题
clear
clc
close all
% 设计参数
f_s = 5000;  % 采样频率 (Hz)
f_c = 800;   % 通带边界频率 (Hz)
f_r = 500;   % 阻带边界频率 (Hz)
delta = 1;   % 通带波动 (dB)
A_t = 40;    % 阻带最小衰减 (dB)
% 使用firpmord估算滤波器阶数
[order, fo, ao, w] = firpmord([f_r f_c], [0 1], [10^(-A_t/20) 1-10^(-delta/20)], f_s);
% 显示估算的滤波器阶数
disp(['估算的滤波器阶数: ', num2str(order)]);

% 使用firpm设计滤波器
b = firpm(order, fo, ao, w);

% 绘制频率响应
[H, f] = freqz(b, 1);  % 计算频率响应

% 绘制幅度响应
figure;
subplot(2,1,1);
plot(f*f_s/2/pi, 20*log10(abs(H)));  % 幅度响应
title('高通滤波器幅度响应');
xlabel('频率 (Hz)');
ylabel('幅度(dB)');
grid on;

% 绘制相位响应
subplot(2,1,2);
plot(f*f_s/2/pi, angle(H));  % 相位响应
title('高通滤波器相位响应');
xlabel('频率 (Hz)');
ylabel('相位 (radians)');
grid on;







%%
% 带宽计算函数
function [bw_3dB, bw_20dB] = calculate_bandwidth(H, F)
    % 计算3dB和20dB带宽
    H_mag = 20 * log10(abs(H));
    
    % 计算3dB带宽
    bw_3dB_start = find(H_mag >= -3, 1, 'first');
    bw_3dB_end = find(H_mag >= -3, 1, 'last');
    bw_3dB = F(bw_3dB_end) - F(bw_3dB_start);

    % 计算20dB带宽
    bw_20dB_start = find(H_mag >= -20, 1, 'first');
    bw_20dB_end = find(H_mag >= -20, 1, 'last');
    bw_20dB = F(bw_20dB_end) - F(bw_20dB_start);
end









