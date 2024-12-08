clc;
clear all;
close all;
%ct
a=[1    -0.886372552937091    0.0752264044191101    0.0518195032914529    0.0356877626508858    0.0245663788107344    0.0168939733876390    0.0115933847488235    0.00792044070864310    0.00535948103002273    0.00355096169035401    0.00224097696680291    0.00396593581833329];%分母
b=0.599222883197722;%分子
Fs = 8000;         % 采样频率 (Hz)
T = 2;             % 信号持续时间 (秒)
N = Fs * T;       % 总采样点数
tt = (0:N-1) / Fs; % 时间向量

% 生成白噪声
rng(1000,'twister');%固定随机数
noise = 10*randn(N, 1);
%ct=lsim(a,b,randn(size(t)));
ct=filter(b,a,noise);%平稳连续谱
ct=ct';
plot(tt,ct);

%lt
% 定义线谱分量的参数
f1 = 97;               % 第一个正弦波频率 (Hz)
A1 = 5;                % 第一个正弦波幅度
phi1=0;
f2 = 36;               % 
f3 = 310;               % 
f5 = 847;               % 
f4 = 1500;               %
% 生成线谱信号
lt = A1 * (sin(2 * pi * f1 * tt + phi1)+sin(2 * pi * f2 * tt + phi1)+sin(2 * pi * f3 * tt + phi1)+sin(2 * pi * f4 * tt + phi1)+sin(2 * pi * f5 * tt + phi1));
plot(tt,lt);


%mt
fp=3;%轴频
nb=5;%叶片数
fb=nb*fp;%叶频
Ai=[0.1,0.1,0.1,0.1,0.16];%调制系数
w0=2*pi*fp;%轴频的角频率
% 
% nn=0:WL-1;
% tt=nn/Fs;
mt=0;
for ii=1:nb
    mt=mt+Ai(ii).*cos(ii*w0.*tt);
end
s=[1+mt].*ct+lt;
s=s';
%plot(mt);
figure;
plot(tt,s);
title("时域图");
xlabel('时间(S)');
ylabel('幅度(V)');
figure;
ffts=fft(s,Fs);
plot(10*log10(abs(ffts.^2)/N/Fs)+120);
xlim([0 Fs/2]);
title("功率谱");
ylabel("幅度/dB");
xlabel("频率/Hz");