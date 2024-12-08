
%% 段仁俊-2024-11-22作业
clear
close all
clc
%% 设置文件夹保存路径
Save_Path='E:\课程\海洋声信号实验\4\Result\';
%% noise
% 生成白噪声
Fs = 5000;         % 采样频率 (Hz)
T = 5;             % 信号持续时间 (秒)
N = Fs * T;       % 总采样点数
tt = (0:N-1) / Fs; % 时间向量
rng(1000,'twister');%固定随机数
noise = randn(N, 1);
%% CW与接收信号Input_signal
%要求：生成一段CW信号。时间长度为1秒，信号到达时间为1.5秒
fs=5000;%采样率
f0=300;%频率
tao=1;%脉宽
t0=1.5;%信号到达时间
t=t0:1/fs:t0+tao-1/fs;
cw=chirp(t,f0,tao,f0);
subplot(2,1,1);
s=[zeros(1.5*Fs,1);cw';zeros(2.5*Fs,1)];
plot(tt,s);
title("CW发射信号时域波形");
subplot(2,1,2);
CwInput_signal=s+noise;
plot(tt,CwInput_signal);
title("接收信号时域波形");
% %% 匹配滤波1
% % 匹配滤波器（脉冲响应为信号的时间反转）
% matched_filter = flipud(s); 
% %plot(tt,matched_filter)
% figure
% % 进行匹配滤波
% Cwoutput_signal = conv(CwInput_signal, matched_filter,'same');
% plot(tt-1,Cwoutput_signal);
%% 匹配滤波2
%相关函数卷积结果
figure
% 使用 xcorr 计算自相关
%[Cwoutput, lags] = xcorr(Input_signal,Input_signal);%如果s与noise不相关可以直接用
[Cwoutput, lags] = xcorr(CwInput_signal,cw');
% 为了绘制输出信号，生成相应的时间轴
lag_time = lags / fs;
plot(lag_time,Cwoutput)
[~,max_time]=max(Cwoutput);
Aim_Length=1500*(max_time-N)/Fs;
title(['目标距离为' num2str(Aim_Length) 'm'])
%% 保存图片
saveas( 1, [Save_Path '第一问CW与发射信号'],'png'); %保存窗口的图像
saveas( 2, [Save_Path '第一问匹配滤波结果'],'png'); %保存窗口的图像

%% LFM信号与接收信号
close all
figure
fs=5000;%采样率
LFMf0=200;%频率
LFMf1=400;%
theta=pi/2;%初始相位
tao=1;%脉宽
LFMt0=2;%信号到达时间
T=5;
N=fs*T;
tt = (0:N-1) / fs; % 时间向量
t=LFMt0:1/fs:LFMt0+tao-1/fs;
LFM=chirp(t,LFMf0,tao,LFMf1);
subplot(2,1,1);
s=[zeros(LFMt0*fs,1);LFM';zeros((T-LFMt0-tao)*fs,1)];
plot(tt,s);
title("LFM发射信号时域波形");
subplot(2,1,2);
LFMInput_signal=s+noise;
plot(tt,LFMInput_signal);
title("接收信号时域波形");
%% 匹配滤波
%相关函数卷积结果
figure
% 使用 xcorr 计算自相关
%[Cwoutput, lags] = xcorr(Input_signal,Input_signal);%如果s与noise不相关可以直接用
[Cwoutput, lags] = xcorr(LFMInput_signal,LFM');
% 为了绘制输出信号，生成相应的时间轴
lag_time = lags / fs;
plot(lag_time,Cwoutput)
[~,max_time]=max(Cwoutput);
Aim_Length=1500*(max_time-N)/Fs;
title(['目标距离为' num2str(Aim_Length) 'm'])
%% 保存图片
saveas( 1, [Save_Path '第二问LFM与发射信号'],'png'); %保存窗口的图像
saveas( 2, [Save_Path '第二问匹配滤波结果'],'png'); %保存窗口的图像



%%
close all
clear 
clc
fs=5000;%采样率
f0=300;%频率
tao=1;%脉宽
t0=1.5;%信号到达时间
t=t0:1/fs:t0+tao-1/fs;
cw=chirp(t,f0,tao,f0);
rng1=cw_noise(5000,5,298,1.5,1);
rng2=cw_noise(5000,5,299,1.5,1);
rng3=cw_noise(5000,5,300,1.5,1);
rng4=cw_noise(5000,5,301,1.5,1);
rng5=cw_noise(5000,5,302,1.5,1);
%% 匹配滤波
%相关函数卷积结果
figure
% 使用 xcorr 计算自相关
[Cwoutput, lags] = xcorr(rng1,cw');
Cwoutput1=xcorr(rng2,cw');
Cwoutput2=xcorr(rng3,cw');
Cwoutput3=xcorr(rng4,cw');
Cwoutput4=xcorr(rng5,cw');
% 为了绘制输出信号，生成相应的时间轴
lag_time = lags / fs;
subplot(2,3,1)
plot(lag_time,Cwoutput)
title('298Hz匹配滤波结果')
subplot(2,3,2)
plot(lag_time,Cwoutput1)
title('299Hz匹配滤波结果')
subplot(2,3,3)
plot(lag_time,Cwoutput2)
title('300Hz匹配滤波结果')
subplot(2,3,4)
plot(lag_time,Cwoutput3)
title('301Hz匹配滤波结果')
subplot(2,3,5)
plot(lag_time,Cwoutput4)
title('302Hz匹配滤波结果')
% [~,max_time]=max(Cwoutput);
% Aim_Length=1500*(max_time-N)/Fs;
% title(['目标距离为' num2str(Aim_Length) 'm'])
%% 保存图片
Save_Path='E:\课程\海洋声信号实验\4\Result\';
saveas( 1, [Save_Path '第三问匹配滤波结果'],'png'); %保存窗口的图像

%%
function CwInput_signal=cw_noise(Fs,T,f0,t0,tao)
% 函数功能为生成CW与noise叠加之后的信号
% Fs 采样频率 (Hz)
% T为噪声总时间 (秒)
%f0为CW信号频率
%t0信号到达时间
%tao为信号持续时间
    N = Fs * T;       % 总采样点数
    rng(1000,'twister');%固定随机数
    noise = randn(N, 1);
    t=t0:1/Fs:t0+tao-1/Fs;
    cw=chirp(t,f0,tao,f0);
    s=[zeros(t0*Fs,1);cw';zeros((T-t0-tao)*Fs,1)];
    CwInput_signal=s+noise;
end




