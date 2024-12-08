%作业
d=0.5;
N=20;
f=5000;
c=1500;
lei=c/f;
theta=-pi/2:0.01:pi/2;
D=sin(N*pi*d*sin(theta)/lei)./(N*sin(pi*d./lei*sin(theta)));
subplot(2,1,1);
plot(theta*180/pi,20*log10(D));
title("d=0.5米");
ylabel('指向性函数/dB');
xlabel('方位角/(°)');
xlim([-90 90]);
ylim([-40 0]);
subplot(2,1,2);
d=0.285;
theta=-pi/2:0.01:pi/2;
D=sin(N*pi*d*sin(theta)/lei)./(N*sin(pi*d./lei*sin(theta)));
plot(theta*180/pi,20*log10(D));
title("d=6/19米");
ylabel('指向性函数/dB');
xlabel('方位角/(°)');
xlim([-90 90]);
ylim([-40 0]);

%%
clc;
clear;
close all;
WL=1024;
STEP=128;
clims=[0,2];
fc=7*10^9;%7GHz
B=4*10^9;%4GHz
Tc=40*10^-6;%40us
S=B/Tc;
fs=2*(fc+B);%采样频率
t=0:1/fs:Tc;
xT=sin(2*pi*(fc*t+0.5*S*t.^2));%chrip信号(TX)
%plot(t,xT);
c=3*10^8;
dmax=fs*c/(2*S);
detad=c/(2*B);
d=10;%墙壁等距离
d1=1;%人距离
delay=2*d/c;
delay1=2*d1/c;
%xR=sin(2*pi*(fc*(t-delay)+0.5*S*(t-delay).^2));%chrip信号(RX)
fif=2*S*d/c;
fif1=2*S*d1/c;
theta=2*pi*fc*delay;
theta1=2*pi*fc*delay1;
fs=10*(fif+fif1);
t=0:1/fs:Tc;
xout=0.3*0.5*cos(2*pi*fif*t+theta)+0.7*0.5*cos(2*pi*fif1*t+theta1);%混频
figure;
fs=2*(fc+B);
spectrogram(xT, hamming(128), 60, 128, fs, 'yaxis');
title('时频域图像TX');
fs=10*fif;
figure;
spectrogram(xout, hamming(128), 60, 128, fs, 'yaxis');
title('有人时频域图像(混频)');
%%
%无人
d1=0;
delay1=2*d/c;
fif1=2*S*d/c;
theta1=2*pi*fc*delay;
xout=0.5*cos(2*pi*fif*t+theta)+0.5*cos(2*pi*fif1*t+theta1);%混频
figure;
spectrogram(xout, hamming(128), 120, 128, fs, 'yaxis');
title('无人时频域图像(混频)');




%[xoutfft,tt,ff]=mySTFT(xout,WL,STEP,fs);
%[xTfft,tt1,ff1]=mySTFT(xT,WL,STEP,fs);
% figure;
% imagesc(tt,ff,abs(xoutfft),clims);  
% axis('xy');
% xlabel('时间 (s)');
% ylabel('频率 (Hz)');
% figure;
% imagesc(tt1,ff1,abs(xTfft),clims);  
% axis('xy');
% xlabel('时间 (s)');
% ylabel('频率 (Hz)');

%%
%2024-10-25
fs=2000;
t1 = 0:1/fs:2;
x1 = chirp(t1,300,2,300,'quadratic');
x2=chirp(t1,300,2,500,'linear');
x3 = chirp(t1,400,2,600,'quadratic');
noise=randn(1,3*length(t1));
x=[x1,x2,x3]+noise;%直接就是一个拼接
nfft=128;
windows=hanning(128);
noverlap=32;%25重叠
spectrogram(x, windows,noverlap, 512,fs);
%%
%
clear;
clc;
close all;
FolderPathAll="E:\海洋声信号实验\ShipRadiatedNoise.mat";
load(FolderPathAll);
fs=5000;
spectrogram(boatdata, hamming(65536/16),64,65536,fs);
%% 2024-10-25作业

%% CW
clear all;clc;
close all
C=1500;
%% 生成信号CW
FS=5e3;
WL=0.5;  
A=1;
T = 0.2;         % 信号脉宽
f0 = 300;       % 信号的初始频率
t = -T/2:1/FS:T/2 ; 
x = A*cos(2*pi*(f0)*t);
%% 

data1=[zeros(1,(WL/2-T/2)*FS) x zeros(1,(WL/2-T/2)*FS)];
%回波接收信号(0,x,0)
% data1=awgn(data1,-5);
% figure(1);plot((1:length(data1))/FS,data1);title('data1');hold on

t=-WL/2:1/FS:WL/2; 
fd=-20:0.2:20; 
for i=1:length(fd) 
    E=exp(-j*2*pi*fd(i).*t);
    Y(i,:)=data1.*E;
end

p=length(data1)*2;
for j=1:length(fd) 
    tmp=fft(xcorr(Y(j,:),data1),p);
    % envelope
    h = [1; 2*ones(fix((p-1)/2),1); ones(1-rem(p,2),1); zeros(fix((p-1)/2),1)];
    M(j,:) = tmp(:).*h;
    M(j,:) = abs(ifft(M(j,:),p))';             
    [Mmax(j)]=max(M(j,:))/p*2;
end
M=M/max(max(M));
t=-WL:1/FS:WL; 
% figure;imagesc(t,fd,M); colorbar;title('信号NBAF'),xlabel('时间 s'),ylabel('频移 Hz');
[tt,ffdd]=meshgrid(t,fd);
figure;mesh(tt,ffdd,M(1:length(fd),1:length(t)));title('信号CW'),xlabel('时间 s'),ylabel('频移 Hz');

%模糊度图及频率时延分辨率
figure;[indxy,indh]=contour(tt,ffdd,M(1:length(fd),1:length(t)),[sqrt(0.5) sqrt(0.5)],'k');title('信号CW-3dB等高线'),xlabel('时间 s'),ylabel('频移 Hz');grid on;% hold on;
h_text=clabel(indxy,indh,'labelSpacing',600);

%% HFM
clear all;clc;
close all
C=1500;
%% 生成信号
FS=5e3;
WL=0.5;  
A=1;
T = 0.2;         % 信号脉宽
t = -T/2:1/FS:T/2; 
tao=0.2;
theta=pi/2;
t0=0;
f1=400;%起始频率
f2=600;%终止频率
k0=(f2-f1)/(tao*f1*f2);
rect_pulse = @(x) double(0<x&x<1);
%chrip中的quadratic好像不是双曲调频，是线性调频
x=A.*rect_pulse((t-t0)/tao).*exp(1i*(-2*pi/k0*log(-k0*t+1/f1)+theta));
%x = A*cos(2*pi*(f0)*t+pi*K*t.^2);
%% 

data1=[zeros(1,(WL/2-T/2)*FS) x zeros(1,(WL/2-T/2)*FS)];
%回波接收信号(0,x,0)
% data1=awgn(data1,-5);
% figure(1);plot((1:length(data1))/FS,data1);title('data1');hold on

t=-WL/2:1/FS:WL/2; 
fd=-20:0.2:20; 
for i=1:length(fd) 
    E=exp(-j*2*pi*fd(i).*t);
    Y(i,:)=data1.*E;
end

p=length(data1)*2;
for j=1:length(fd) 
    tmp=fft(xcorr(Y(j,:),data1),p);
    % envelope
    h = [1; 2*ones(fix((p-1)/2),1); ones(1-rem(p,2),1); zeros(fix((p-1)/2),1)];
    M(j,:) = tmp(:).*h;
    M(j,:) = abs(ifft(M(j,:),p))';             
    [Mmax(j)]=max(M(j,:))/p*2;
end
M=M/max(max(M));
t=-WL:1/FS:WL; 
% figure;imagesc(t,fd,M); colorbar;title('信号NBAF'),xlabel('时间 s'),ylabel('频移 Hz');
[tt,ffdd]=meshgrid(t,fd);
figure;mesh(tt,ffdd,M(1:length(fd),1:length(t)));title('信号HFM'),xlabel('时间 s'),ylabel('频移 Hz');

%模糊度图及频率时延分辨率
figure;[indxy,indh]=contour(tt,ffdd,M(1:length(fd),1:length(t)),[sqrt(0.5) sqrt(0.5)],'k');title('信号HFM-3dB等高线'),xlabel('时间 s'),ylabel('频移 Hz');grid on;% hold on;
h_text=clabel(indxy,indh,'labelSpacing',600);

%% 组合
clear all;clc;
close all
C=1500;
%% 生成信号
FS=5e3;
WL=1;  
A=1;
T = 0.4;         % 信号脉宽
f0 = 300;       % CW信号的初始频率
t1 =-T/2:1/FS:0 ; 
t2 =1/FS:1/FS:T/2 ;
t=-T/2:1/FS:T/2;
tao=0.2;
theta=pi/2;
t0=0;
f1=400;%起始频率
f2=600;%终止频率
k0=(f2-f1)/(tao*f1*f2);
rect_pulse = @(x) double(0<x&x<1);
%chrip中的quadratic好像不是双曲调频，是线性调频
x1=A.*rect_pulse((t1-t0)/tao).*exp(1i*(-2*pi/k0*log(-k0*t1+1/f1)+theta));
x2 = A*cos(2*pi*(f0)*t2);
x=[x1,x2];
%% 

data1=[zeros(1,(WL/2-T/2)*FS) x zeros(1,(WL/2-T/2)*FS)];
%回波接收信号(0,x,0)
% data1=awgn(data1,-5);
% figure(1);plot((1:length(data1))/FS,data1);title('data1');hold on

t=-WL/2:1/FS:WL/2; 
fd=-20:0.4:20; 
for i=1:length(fd) 
    E=exp(-j*2*pi*fd(i).*t);
    Y(i,:)=data1.*E;
end

p=length(data1)*2;
for j=1:length(fd) 
    tmp=fft(xcorr(Y(j,:),data1),p);
    % envelope
    h = [1; 2*ones(fix((p-1)/2),1); ones(1-rem(p,2),1); zeros(fix((p-1)/2),1)];
    M(j,:) = tmp(:).*h;
    M(j,:) = abs(ifft(M(j,:),p))';             
    [Mmax(j)]=max(M(j,:))/p*2;
end
M=M/max(max(M));
t=-WL:1/FS:WL; 
% figure;imagesc(t,fd,M); colorbar;title('信号NBAF'),xlabel('时间 s'),ylabel('频移 Hz');
[tt,ffdd]=meshgrid(t,fd);
figure;mesh(tt,ffdd,M(1:length(fd),1:length(t)));title('组合信号NBAF'),xlabel('时间 s'),ylabel('频移 Hz');

%模糊度图及频率时延分辨率
figure;[indxy,indh]=contour(tt,ffdd,M(1:length(fd),1:length(t)),[sqrt(0.5) sqrt(0.5)],'k');title('组合信号-3dB等高线'),xlabel('时间 s'),ylabel('频移 Hz');grid on;% hold on;
h_text=clabel(indxy,indh,'labelSpacing',600);








%%
function [yf,tt,ff]=mySTFT(noisy_sig,wLen,step,fs)
% 计算含噪信号的分帧功率谱/周期图
overlap_factor=max(wLen/step,1);            % 重叠系数=窗长/步进
t_increment=step/fs;                        % 所需帧增量
n_increment=pow2(nextpow2(step*sqrt(0.5))); % 将步进step转化为2的幂次
fftLen=n_increment*overlap_factor;          % FFT长度
w=sqrt(hamming(fftLen+1))';                 % 使用sqrt hamming window
w(end)=[];
w=w/sqrt(sum(w(1:n_increment:fftLen).^2));  % 归一化 normalize to give overall gain of 1

% 分帧 开始
noisyLen=length(noisy_sig);                                           % 含噪信号的长度
n_wLen=length(w);                                                     % 窗长
w=w(:).';
n=noisyLen-n_wLen+n_increment;                                        % 进行分帧的点数=总数据长度-窗长+帧增量
num_frame=max(fix(n/n_increment),0);                                  % 总帧数
num_left=n-n_increment*num_frame+(num_frame==0)*(n_wLen-n_increment); % 分完剩余的采样点数
if num_left>0
    left_flag=1;                                                      % 如果有剩余的采样点，再补一帧
else
    left_flag=0;
end
f=zeros(num_frame+left_flag,n_wLen);
indf=n_increment*(0:(num_frame-1)).';
inds=(1:n_wLen);
if left_flag==1                                                       % 类似于将数据长度凑成可以正好分帧的帧数
    f(1:num_frame,:)=noisy_sig(indf(:,ones(1,n_wLen))+inds(ones(num_frame,1),:));
    ix=1+mod((num_frame*n_increment):(num_frame*n_increment+n_wLen-1),2*noisyLen);
    f(num_frame+1,:)=noisy_sig(ix+(ix>noisyLen).*(2*noisyLen+1-2*ix));
    num_frame=size(f,1);
else
    f(:)=noisy_sig(indf(:,ones(1,n_wLen))+inds(ones(num_frame,1),:)); % 
end
f=f.*w(ones(num_frame,1),:);
t0=(1+n_wLen)/2;
t=t0+n_increment*(0:(num_frame-1)).';
tt=t/fs; % 时间轴
% 分帧 结束

% 计算分帧fft 开始
s=size(f);
if prod(s)==1
    yf=f; % 如果只有一个数据，不用做FFT，yf为含噪信号周期图
else
    d=2;
    yf=fft(f,fftLen,d);
    yf=reshape(yf,prod(s(1:d-1)),fftLen,prod(s(d+1:end)));
    s(d)=1+fix(fftLen/2);
    yf(:,s(d)+1:end,:)=[];
    yf=reshape(yf,s);
end
% 计算分帧fft 结束

[nr,nf2]=size(yf);
yf(1,1,1)=yf(1,1);
yf=reshape(yf,nr,nf2);

% !!!

% yf=newyf;
yp=yf(:,1:nf2).*conj(yf(:,1:nf2));
ff=(0:nf2-1)*fs/fftLen; % 频率轴
end











