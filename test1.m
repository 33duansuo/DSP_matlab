clc;
clear;
fs=8000;%采样率
N=fs;
A=1;%仿真幅度
f0=500;%中心频率
theta=pi/2;%初始相位
tao=1;%脉宽
t0=0;%信号到达时间
rect_pulse = @(x) double(0<x&x<1);
t=0:1/fs:1;
s=A.*rect_pulse((t-t0)/tao).*exp(1i*(2*pi*f0*t+theta));
ffts=fft(s,N);
ft=(0:N-1)*fs/N;
subplot(2,1,1);
plot(t,s);
title("CW时域");
subplot(2,1,2);
plot(ft,abs(ffts));
title("CW频域");
%%
%LFM
clc;
clear;
close all;
A=1;
tao=1;
theta=pi/2;
t0=0;
fs=8000;
f1=200;
u=300;
t=0:1/fs:1;
N=fs;
rect_pulse = @(x) double(0<x&x<1);
s=A.*rect_pulse((t-t0)/tao).*exp(1i*(2*pi*(0.5*u*t.^2+f1*t)+theta));
ffts=fft(s,N);
ft=(0:N-1)*fs/N;
subplot(2,1,1);
plot(t,s);
title("LFM时域");
subplot(2,1,2);
plot(ft,abs(ffts));
xlim([200 300]);
title("LFM频域");
%%
%HFM
clc;
clear;
close all;
A=1;
tao=1;
theta=pi/2;
t0=0;
fs=8000;
f1=200;%起始频率
f2=500;%终止频率
k0=(f2-f1)/(tao*f1*f2);
t=0:1/fs:1;
N=fs;
rect_pulse = @(x) double(0<x&x<1);
s=A.*rect_pulse((t-t0)/tao).*exp(1i*(-2*pi/k0*log(-k0*t+1/f1)+theta));
ffts=fft(s,N);
ft=(0:N-1)*fs/N;
subplot(2,1,1);
plot(t,s);
title("HFM时域");
subplot(2,1,2);
plot(ft,abs(ffts));
title("HFM频域");
xlim([200 500]);

%%
%4
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
%% 阵列延时
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
%ffts=fft(snew,Fs);
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
%% 波束形成
close all;
M=16;
d=1;
c=1500;
filterOrder=10;
zeroslength=(M-1)*d/c*Fs;
theta1=(0:1:180)*pi/180;%扫描角度
BeamNum=length(theta1);
BeamEnergy=zeros(1,BeamNum);
WL=length(xxdelay)+2*zeroslength;
 for ii=1:BeamNum%波束
     Yp=zeros(1,WL);%预成波束数据
    for jj=1:M%阵元
        DFS=-(jj-1)*d*cos(theta1(ii))/c*Fs;
        DI=DFS-round(DFS);
        outData=Delayfilter(xxdelay(jj,:),DFS,filterOrder,zeroslength);
        Yp=Yp+outData;%确保是从0时刻开始
    end
    Yp=Yp/M;
    % [b,a]=butter(10,[350/(Fs/2) 750/(Fs/2)],'bandpass');
    % Yp11=filter(b,a,Yp1);
     YpFFT=fft(Yp)/WL;
    BeamEnergy(ii)=sum(abs(YpFFT).^2,'all');
 end  
figure(2)
plot(theta1/pi*180,10*log10(BeamEnergy/max(BeamEnergy)),'b-','LineWidth',1)
xlabel('\fontsize{15}angle (°)')
ylabel('\fontsize{15}Power (dB)')
title('60°入射信号波束扫描结果')
















% %%
% 
% function [b, a] = frac_delay_filter(N, d)
%     % 创建窗口函数（Hamming窗）
%     w = hamming(N+1);
% 
%     % 计算小数时延滤波器的冲激响应
%     b = zeros(1, N+1);
%     for n = 0:N
%         b(n+1) = sin(pi * (n-N/2-d)) / (pi * (n-N/2-d));
%         b(n+1) = b(n+1)*w(n+1); % 乘以窗口
%     end
% 
%     % 归一化使得总增益为1
%     b = b / sum(b);
%     a = 1; % FIR滤波器，所以a=1
% end



