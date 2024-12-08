clc;
clear;
close all;
pic1=figure;
pic2=figure;
[y,Fs]=audioread("D:\BaiduNetdiskDownload\2024第1周实验要求与数据\data\R1.wav");%32000Hz
t=0:1/32000:(960000-1)*1/32000;%数据时长
fprintf('舰船辐射噪声数据时长(秒): %.2f\n', 1/Fs*length(y)); % 保留两位小数
fprintf('舰船辐射噪声数据点数: %.2f\n', length(y)); 
fprintf('舰船辐射噪声采样率(Hz): %.2f\n', Fs); 
[y1,Fs1]=audioread("D:\BaiduNetdiskDownload\2024第1周实验要求与数据\data\P1.wav");%10000Hz
t1=0:1/10000:(50000-1)*1/10000;%数据时长
fprintf('脉冲信号数据时长(秒): %.2f\n', 1/Fs1*length(y1)); % 保留两位小数
fprintf('脉冲信号数据点数: %.2f\n', length(y1)); 
fprintf('脉冲信号采样率(Hz): %.2f\n', Fs1); 
figure(pic1);
subplot(2,1,1);
plot(t,y,'r');
title("舰船辐射噪声时域")
subplot(2,1,2);
plot(t1,y1,'b');
title("脉冲时域")
%%
%播放
sound(y,Fs);
sound(y1,Fs1);
%%
%统计

M=mean(y);
fprintf('舰船辐射噪声均值: %.2f\n', M); % 保留两位小数
M1=mean(y1);
fprintf('脉冲信号均值: %.2f\n', M1);
S=std(y);
fprintf('舰船辐射噪声标准差: %.2f\n', S); % 保留两位小数
S1=std(y1);
fprintf('脉冲信号标准差: %.2f\n', S1); % 保留两位小数
V=var(y);
fprintf('舰船辐射噪声方差: %.2f\n', V); % 保留两位小数
V1=var(y1);
fprintf('脉冲信号方差: %.2f\n', V1); % 保留两位小数

%%
%功率谱
nfft=Fs/10;%默认输入信号的长度
nfft1=Fs1/10;
figure(pic2);
[Pxx,F]=periodogram(y,[],nfft,Fs);%舰船辐射噪声信号
[Pxx1,F1]=periodogram(y1,[],nfft1,Fs1);%脉冲信号
subplot(2,1,1);
Pxx(1)=Pxx(2);%对0进行操作
Pxx1(1)=Pxx1(2);
plot(F,10*log10(Pxx));%功率谱
title("舰船辐射噪声功率谱")
subplot(2,1,2);
plot(F1,10*log10(Pxx1));
title("脉冲功率谱")
%%
%时频图
[s,f,t]=spectrogram(y,256,128,2048,Fs,'yaxis');
[s1,f1,t1]=spectrogram(y1,256,128,2048,Fs1,'yaxis');
figure;
imagesc(t,f,10*log10(abs(s)));
title("舰船辐射噪声时频图")
xlabel('Time (s)');
ylabel('Frequency (Hz)');
figure;
imagesc(t1,f1,10*log10(abs(s1)));
title("脉冲时频图")
xlabel('Time (s)');
ylabel('Frequency (Hz)');