clc;
clear;
close all;
% 参数设置
N = 10; % 滤波器阶数
delay = 0.5; % 小数延时

[b, a] = frac_delay_filter(N, delay);

% 计算频率响应
[H, f] = freqz(b,a); 

% 绘制幅频响应和相频响应
figure;

% 幅频响应
subplot(2, 1, 1);
plot(f/pi, (abs(H)));
title("幅度响应”");
xlabel("归一化频率");
ylabel("幅度（dB）");
grid on;

% 相频响应s
subplot(2, 1, 2);
plot(f/pi, unwrap(angle(H)*180/pi));
title("相位响应");
xlabel("归一化频率");
ylabel('相位(°)');
grid on;



figure;
f=500;
fs=4000;
t=0:N;
x=sin(2*pi*t*f/fs);
x1=sin(2*pi*(t-0.5)*f/fs);
y=conv(x,b,'same');
plot(t,y,'r');
hold on;
plot(t,x,'b','Marker','o');
hold on;
plot(t,x1,'g','Marker','+');
grid on;
xlabel("时间序列索引");
ylabel('信号');
legend("滤波器","原信号","实际延时信号");
%%
function [h,a]=frac_delay_filter(N, tao)
M=N/2;

 k=-M:M;
    hd=sinc(k-tao);
    hnWin=hanning(2*M+1)';
    hd=hd.*hnWin;
    h=hd/sum(hd);
    a=1;

end

