close all
clear all
clc
%%


% 加载数据
Data = load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\TraceBeamData1.mat');
TraceBeamData = Data.TraceBeamData1;
%%

% 参数设置
Fs = 20000;           % 采样率 20 kHz
WL = 2 * Fs;          % 每个窗口的长度
STEP = 0.5 * WL;      % 步长，50% 重叠
Nfft = WL;            % FFT点数与窗口长度相等
Noverlap = STEP;      % 重叠长度
TraceBeamPsd1 = zeros(1, WL);  % 初始化功率谱

% 带通滤波器设计
startFreq = 50;        % 分析频带频率下限
endFreq = 820;         % 分析频带频率上限
bpFilterOrder = 4096;  % 带通滤波器阶数
wp = [startFreq / (Fs / 2), endFreq / (Fs / 2)];  % 归一化截止频率
b = fir1(bpFilterOrder, wp);  % 设计带通滤波器

% Welch方法计算功率谱并应用带通滤波器
for ii = 1:16
    data = TraceBeamData((ii-1)*STEP + (1:WL));  % 按照步长取数据
    
    % 应用带通滤波器
    filteredData = filter(b, 1, data);  % 使用带通滤波器
    
    % 计算滤波后数据的功率谱
    Pxx = abs(fft(filteredData, Nfft)).^2 / WL;
    TraceBeamPsd1 = TraceBeamPsd1 + Pxx;  % 累加功率谱
end

% 求平均功率谱
TraceBeamPsd1 = TraceBeamPsd1 / 16;  % 取平均值
%% 线谱检测
% 计算频率轴
fPsd = (0:WL-1) * Fs / Nfft;  % 频率分辨率

% 连续谱平滑（使用中值滤波器）
smoothedPsd = medfilt1(TraceBeamPsd1, 50);  % 50点中值滤波

% 计算差值谱（即背景均衡后的功率谱）
residualPsd = TraceBeamPsd1 - smoothedPsd;  % 差值谱

% 线谱检测（使用 3σ 标准来识别线谱）
meanRes = mean(residualPsd);   % 差值谱均值
stdRes = std(residualPsd);     % 差值谱标准差
threshold = meanRes + 3*stdRes.^2;  % 线谱检测阈值
cnt=1;
% 找出所有超过阈值的点
lineSpectra = find(residualPsd > threshold);
for i=2:length(lineSpectra)-1
    if residualPsd(lineSpectra(i))<residualPsd(lineSpectra(i-1))||residualPsd(lineSpectra(i))<residualPsd(lineSpectra(i+1))
        delete_index(cnt)=i;
        cnt=cnt+1;
    end
end
lineSpectra(delete_index)=[];
% 绘制结果
figure(2);
plot(fPsd, 10*log10(TraceBeamPsd1));  % 绘制原始功率谱
axis([50 820 0 50])
hold on;
plot(fPsd(lineSpectra), 10*log10(TraceBeamPsd1(lineSpectra)), 'ro');  % 标记线谱点
xlabel('频率(Hz)');
ylabel('功率谱密度 (dB/Hz)');
title('线谱检测');
legend('功率谱密度', '线谱检测');

% 绘制背景均衡后的功率谱（差值谱）
figure(3);
plot(fPsd, 10*log10(residualPsd));  % 绘制差值谱
xlabel('Frequency (Hz)');
ylabel('Residual Power (dB)');
title('差值谱');
axis([50 820 0 50])
%%
close all
clear all
clc
% 加载数据
Data = load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\TraceBeamData1.mat');
TraceBeamData = Data.TraceBeamData1;

%% 解调谱
% 参数设置
Fs = 20000;                     % 采样率
WL = 2 * Fs;          % 每个窗口的长度
STEP = 0.5 * WL;      % 步长，50% 重叠
Noverlap = STEP;      % 重叠长度
TraceBeamPsd1 = zeros(1, WL);  % 初始化功率谱
freqBands =linspace(1000,5000,10);
lowCutoffFreq = 50;              % 低通滤波器截止频率
lpFilterOrder = 10;             % 低通滤波器阶数
Nfft = 2 * Fs;                   % FFT 点数
Fmin = 1;                        % 轴频搜索最小值
Fmax = 20;                       % 轴频搜索最大值
L = 10;                          % 搜索的倍频次

% 初始化结果变量
subBandSpectra = zeros(length(freqBands)-1, 10*Nfft);  % 存储每个子频带的解调谱

% 逐个频带进行分析
for i = 1:length(freqBands)-1
    % 当前频带
    bandStart = freqBands(i);
    bandEnd = freqBands(i+1);
    % 带通滤波器设计
    wp = [bandStart / (Fs / 2), bandEnd / (Fs / 2)];  % 归一化带通频率
    b_bp = fir1(4096, wp);                           % 带通滤波器设计
    % for ii = 1:16
    % data = TraceBeamData((ii-1)*STEP + (1:WL));  % 按照步长取数据
    % % 应用带通滤波器
    % filteredData = filter(b_bp, 1, data);  % 使用带通滤波器
    % TraceBeamPsd1=filteredData+TraceBeamPsd1;
    % end
    % 平方运算
    TraceBeamPsd1 = filter(b_bp, 1, TraceBeamData);  % 使用带通滤波器
    squaredSignal = TraceBeamPsd1.^2;               % 对带通滤波后的信号进行平方
    % 低通滤波器设计
    wl =lowCutoffFreq / (Fs / 2);                   % 归一化低通截止频率
    b_lp = fir1(lpFilterOrder, wl);                  % 低通滤波器设计

    % 对平方信号进行低通滤波
    envelopeSignal=filter(b_lp, 1, squaredSignal); % 提取包络信号
    Pxx = abs(fft(envelopeSignal, Nfft*10)).^2 /WL;  % 计算功率谱

    % 对包络信号进行谱估计
    subBandSpectra(i, :) = Pxx;                     % 保存当前频带的解调谱
end

% 分频带融合
fusedSpectrum = sum(subBandSpectra, 1);             % 融合所有子频带的解调谱
f = (0:10*Nfft-1) * Fs / (10*Nfft);                         % 频率轴

% 绘制融合后的解调谱
figure;
plot(f, 10*log10(fusedSpectrum));
xlabel('Frequency (Hz)');
ylabel('PSD(dB/Hz)');
title('解调谱');
axis([0 50 0 40])
% % 轴叶频提取
% % 限制频率范围到 [Fmin, Fmax]
% freqSearchRange = (f >= Fmin) & (f <= Fmax);
% searchFrequencies = f(freqSearchRange);
% 
% % 轴频检测
% D = zeros(size(searchFrequencies));  % 初始化轴频指标
% for j = 1:length(searchFrequencies)
%     fj = searchFrequencies(j);
%     harmonicIndices = round(fj * (1:L) * 10*Nfft / Fs);  % 倍频对应的频率索引
%     harmonicIndices = harmonicIndices(harmonicIndices <= 10*Nfft);  % 剔除超出范围的频率
%     D(j) = sum(fusedSpectrum(harmonicIndices)) / L;   % 计算轴频指标
% end
% 
% % 找到最大 D 对应的轴频
% [~, maxIndex] = max(D);
% shaftFrequency = searchFrequencies(maxIndex);  % 最大 D 对应的轴频
% 
% % 绘制轴频搜索结果
% figure;
% plot(searchFrequencies, D, '-o');
% xlabel('Frequency (Hz)');
% ylabel('D_j');
% title('Shaft Frequency Detection');
% grid on;
% % 输出结果
% disp(['Detected Shaft Frequency: ', num2str(shaftFrequency), ' Hz']);












