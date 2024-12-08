
%% 舰船辐射噪声
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

%% 阵列延时
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















