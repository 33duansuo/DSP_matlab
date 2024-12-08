clear
clc
close all
%% 程序版本信息
% 作者： 段仁俊
% 时间： 20241130
%% 设定数据文件保存地址
FolderPath='E:\课程\DSP\';
% 数据文件名
Save_Path = strcat(FolderPath,'Result1\');
% 处理结果保存地址
%% ===================【程序说明】===================
% 无
%% 第一题(模拟滤波器)
close all
fc=300;%通带边界频率
fr=200;%阻带边界频率
fs=1000;%采样频率
Rp=0.8;%dB
Rs=20;%dB
Wp=fc*2*pi;%通带角频率
Ws=fr*2*pi;%阻带角频率
[N,Wn]=cheb1ord(Wp,Ws,Rp,Rs,'s');
[B,A]=cheby1(N,Rp,Wn,'high','s');
omega=[Wp Ws];
h1 = freqs(B,A,omega);
Ap=-20*log10(abs(h1(1)));
As=-20*log10(abs(h1(2)));
[h,w]=freqs(B,A);
f=w/(2*pi);
plot(f,20*log(abs(h)))
axis([100 800 -200 0])
xlabel('频率/Hz'); ylabel('幅度/dB');
title(['切比雪夫高通滤波器 Ap=' num2str(Ap) 'dB' ' As=' num2str(As) 'dB'])
%% 保存图片
saveas(1, [Save_Path '第一问'],'svg'); 
close all;
%% 第二题
close all
clear
clc
fc=200;%通带边界频率
fr=300;%阻带边界频率
fs=1000;%采样频率
Rp=1;%dB
Rs=25;%dB

% 脉冲响应不变法
Wp=fc*2*pi;%通带角频率
Ws=fr*2*pi;%阻带角频率
[N,Wn]=buttord(Wp,Ws,Rp,Rs,'s');
[B,A]=butter(N,Wn,'low','s');
[BZ,AZ]=impinvar(B,A,fs);

% 双线性预畸变
Wp1=2*fs*tan(2*pi*fc/(2*fs));
Ws1=2*fs*tan(2*pi*fr/(2*fs));
[N1,Wn1]=buttord(Wp1,Ws1,Rp,Rs,'s');
[B1,A1]=butter(N1,Wn1,'low','s');
[BZ1,AZ1]=bilinear(B1,A1,fs);%双线性

%记录衰减量和带宽
omega=[Wp/(fs) Ws/(fs)];
h11 = freqz(BZ,AZ,omega);
Ap=-20*log10(abs(h11(1)));
As=-20*log10(abs(h11(2)));%脉冲响应不变衰减量

h12 = freqz(BZ1,AZ1,omega);
Ap1=-20*log10(abs(h12(1)));
As1=-20*log10(abs(h12(2)));%双线性衰减量

[h1,w]=freqz(BZ,AZ);
[h2,w1]=freqz(BZ1,AZ1);
f=w/(2*pi)*fs;
f1=w1/(2*pi)*fs;
plot(f,20*log(abs(h1)),'-.',f1,20*log(abs(h2)),'-','LineWidth',1.5)
xlabel('频率/Hz'); ylabel('幅度/dB');
legend('脉冲响应不变','双线性')
title('巴特沃斯数字低通滤波器频率响应')
y1=20*log(abs(h1));
y2=20*log(abs(h2));
[~,dBbandwith]=min(abs(y1(1:end)+3));
[~,dBbandwith1]=min(abs(y2(1:end)+3));
disp(['【脉冲响应不变法的3dB带宽为】' num2str(dBbandwith) 'Hz' ' 衰减量为' num2str(As) 'dB'])
disp(['【双线性变换法的3dB带宽为】' num2str(dBbandwith1) 'Hz' ' 衰减量为' num2str(As1) 'dB'])
%% 保存图片
saveas(1, [Save_Path '第二问'],'svg'); 
close all;
%% 第三问
close all
clear all
clc
fc=1200;%通带边界频率
fr=2000;%阻带边界频率
fs=8000;%采样频率
Rp=0.5;%dB
Rs=40;%dB

% 双线性预畸变(巴特沃斯型)
Wp=2*fs*tan(2*pi*fc/(2*fs));
Ws=2*fs*tan(2*pi*fr/(2*fs));
[N,Wn]=buttord(Wp,Ws,Rp,Rs,'s');
[B,A]=butter(N,Wn,'low','s');
[BZ,AZ]=bilinear(B,A,fs);%双线性

% 双线性预畸变(切比雪夫型)
[N1,Wn1]=cheb1ord(Wp,Ws,Rp,Rs,'s');
[B1,A1]=cheby1(N1,Rp,Wn1,'low','s');
[BZ1,AZ1]=bilinear(B1,A1,fs);%双线性

% 双线性预畸变(椭圆型)
[N2,Wn2]=ellipord(Wp,Ws,Rp,Rs,'s');
[B2,A2]=ellip(N2,Rp,Rs,Wn2,'low','s');
[BZ2,AZ2]=bilinear(B2,A2,fs);%双线性

[h1,w]=freqz(BZ,AZ);
[h2,~]=freqz(BZ1,AZ1);
[h3,~]=freqz(BZ2,AZ2);
f=w/(2*pi)*fs;
figure
plot(f,20*log(abs(h1)))
title(['巴特沃斯型 阶数=' num2str(N)])
xlabel('频率/Hz'); ylabel('幅度/dB');
figure
plot(f,20*log(abs(h2)))
xlabel('频率/Hz'); ylabel('幅度/dB');
title(['切比雪夫型 阶数=' num2str(N1)])
figure
plot(f,20*log(abs(h3)))
xlabel('频率/Hz'); ylabel('幅度/dB');
title(['椭圆型 阶数=' num2str(N2)])
%% 保存图片
saveas(1, [Save_Path '第三问0'],'svg'); 
saveas(2, [Save_Path '第三问1'],'svg'); 
saveas(3, [Save_Path '第三问2'],'svg'); 
close all;
%% 第四问
close all
clear all
clc
fs=30000;
f1=2000;
f2=3000;
f3=1500;%下阻带
f4=6000;%上阻带
Rp=5;
Rs=20;
%双线性法
w1=2*fs*tan(2*pi*f1/(2*fs));%预畸变
w2=2*fs*tan(2*pi*f2/(2*fs));
w3=2*fs*tan(2*pi*f3/(2*fs));
w4=2*fs*tan(2*pi*f4/(2*fs));
[n, Wn] = buttord([w1 w2],[w3 w4],3,20,'s'); 
[B,A] = butter(n, Wn, 's');
[BZ,AZ]=bilinear(B,A,fs);
[h,w]=freqz(BZ,AZ);
f=w/(2*pi)*fs;
%脉冲响应不变法
[n1, Wn1] = buttord([f1 f2]*2*pi,[f3 f4]*2*pi,3,20,'s'); 
[B1,A1] = butter(n1, Wn1, 's');
[BZ1,AZ1]=impinvar(B1,A1,fs);
[h1,~]=freqz(BZ1,AZ1);
plot(f,20*log10(abs(h)),'-.',f,20*log10(abs(h1)),'-','LineWidth',1.5);
grid; xlabel('频率/Hz'); ylabel('幅度');
legend('双线性法','脉冲响应不变法')
axis([1000 7000 -100 0])
title('巴特沃斯带通幅频响应')
%% 保存图片
saveas(1, [Save_Path '第四问'],'svg'); 
close all;
%% 第五问
close all
clear all
clc
fs=10000;
f1=1000;
f2=2000;
f3=500;%下通带
f4=3000;%上通带
%双线性法
Rp=3;
Rs=18;
w1=2*fs*tan(2*pi*f1/(2*fs));%预畸变
w2=2*fs*tan(2*pi*f2/(2*fs));
w3=2*fs*tan(2*pi*f3/(2*fs));
w4=2*fs*tan(2*pi*f4/(2*fs));
[n, Wn] = cheb1ord([w3 w4],[w1 w2],Rp,Rs,'s'); 
[B,A] = cheby1(n,Rp,Wn,'stop','s');
[BZ,AZ]=bilinear(B,A,fs);
[h,w]=freqz(BZ,AZ);
f=w/(2*pi)*fs;
plot(f,20*log10(abs(h)),'LineWidth',1.5)
grid; xlabel('频率/Hz'); ylabel('幅度');
title('切比雪夫带阻幅频响应')
%% 保存图片
saveas(1, [Save_Path '第五问'],'svg'); 
close all;
