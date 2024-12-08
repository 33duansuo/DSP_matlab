%% 阵列延时
clc;
close all;
clear all;
Fs=8000;
T=2;
N=Fs*T;
tt = (0:N-1) / Fs; % 时间向量
rng(1000,"twister");
s=randn(N,1);
close all;
d=1;
theta0=60/180*pi;%目标方位
c=1500;%水中声速
MM=16;%阵元个数
for m=1:MM
tm=(m-1)*d*cos(theta0)/c*Fs;%数字整体时延
outdata=Delayfilter(s',tm,10,70);%70为补0的个数
xxdelay(m,:)=outdata;
end
TT=(0:length(xxdelay)-1)/Fs;%更新后的时间向量
%% 作图
close all;
figure;
plot(tt,s,'r');
hold on;
plot(TT,xxdelay(16,:),'b');
legend("原舰船辐射噪声","线列阵接收后信号");
title("线列阵接收后时域图");
xlabel("时间/s");
ylabel("幅度");
figure;
snew=s;
ffts=fft(snew,Fs);
plot(10*log10(abs(ffts.^2)/N/Fs)+120,'r');
hold  on;
xlim([0 Fs/2]);
title("功率谱");
ylabel("幅度/dB");
xlabel("频率/Hz");

fft1=fft(xxdelay(16,:),2*Fs);
xlim([0 Fs/2]);
plot(10*log10(abs(fft1.^2)/N/Fs)+120,'b');
title("线列阵接收后功率谱");
legend("原舰船辐射噪声","线列阵接收后信号");
ylabel("幅度/dB");
xlabel("频率/Hz");
%% 延时图
figure
imagesc(tt,1:m,xxdelay);