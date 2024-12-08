clear
clc
close all
%% 程序版本信息
% 作者： 段仁俊
% 时间： 20241114
%% 设定数据文件保存地址
FolderPath='E:\课程\DSP\';
% 数据文件名
Save_Path = strcat(FolderPath,'Result\');
% 处理结果保存地址
%% ===================【程序说明】===================
%参考文献[1]张鹏举,戚晨皓.基于线性相关的不同点循环相关计算方法[J].电气电子教学学报,2022,44(06):129-132.
%参考文献[2]吴镇扬．数字信号处理[M].第3版.北京:高等教育出版社,2016．
%若需保存图片，合理安排文件夹路径
%程序请按照分节的顺序运行，无需额外清空变量
%  调用子程序：  
%  （1）overlap_add_convolution重叠相加法
%  （2）overlap_keep_convolution重叠保留法
%% 高斯序列
Fs=1;
n=0:1/Fs:15;
p=[8 13 14];
q=[2 4 8];
%固定p,改变q
for i=1:length(q)
   x_a(i,:)=exp(-(n-p(1)).^2./q(i));
end
for i=1:length(q)
    subplot(3,1,i);
    stem(n,x_a(i,:));
    if i==length(q)%美观
        xlabel("时间索引")
        ylabel("幅度")
    end
    title(['时域p=8,q=' num2str(2+2*i)])

end
figure;
for j=1:length(q)
    subplot(3,1,j);
    x_afft(j,:)=real(fft(x_a(j,:)));
    stem(abs(x_afft(j,:)));
    if j==length(q)%美观
        xlabel("频率索引")
        ylabel("幅度")
    end
    title(['频域p=8,q=' num2str(2+2*j)])
end
%固定q,改变p
for i=1:length(p)
   x_a(i,:)=exp(-(n-p(i)).^2./q(3));
end
figure;
for i=1:length(p)
    subplot(3,1,i);
    stem(n,x_a(i,:));
    if i==length(q)%美观
        xlabel("时间索引")
        ylabel("幅度")
    end
    title(['时域q=8,p=' num2str(p(i))])

end
figure;
for j=1:length(q)
    subplot(3,1,j);
    x_afft(j,:)=real(fft(x_a(j,:)));
    stem(abs(x_afft(j,:)));
    if j==length(q)%美观
        xlabel("频率索引")
        ylabel("幅度")
    end
    title(['频域q=8,p=' num2str(p(j))])
end



%% 保存图片
%需要设置保存路径，不保存图片无需运行
saveas( 1, [Save_Path '时域P不变'],'svg'); 
saveas( 2, [Save_Path '频域P不变'],'svg');
saveas( 3, [Save_Path '时域q不变'],'svg'); 
saveas( 4, [Save_Path '频域q不变'],'svg'); 
close all;
%% 第二问
close all
clc
%% 衰减正弦序列
Fs=1;
n=0:1/Fs:15;
a=0.1;
f=[0.0625 0.4375 0.5625];
for i=1:length(f)
    x_b(i,:)=exp(-a*n).*sin(2*pi*f(i)*n);
    subplot(3,2,2*(i-1)+1);
    stem(n,x_b(i,:));
    xlabel("时间索引")
    ylabel("幅度")
    title(['a=0.1,f=' num2str(f(i))])
end
for j=1:length(f)
    subplot(length(f),2,2*j);
    x_bfft(j,:)=fft(x_b(j,:));
    stem(n,abs(x_bfft(j,:)));
    xlabel("频率索引")
    ylabel("幅度")
    title(['a=0.1' 'f=' num2str(f(j))])
end
%% 保存图片
saveas( 1, [Save_Path '衰减正弦f变化'],'svg'); 
close all;
%% 第三问三角波和反三角波
close all
for n=0:7
    if n<=4
        x_c(n+1)=n;
    end
    if n>4
        x_c(n+1)=8-n;
    end
end
figure;
subplot(2,3,1)
stem(0:length(x_c)-1,x_c);
axis([0 8 0 5]);
title('三角波')
subplot(2,3,2)
x_cfft=fft(x_c,8);
stem(0:length(x_cfft)-1,abs(x_cfft));
axis([0 8 0 20]);
title('三角波8点FFT')
subplot(2,3,3)
x_c1fft=fft(x_c,32);
stem(0:length(x_c1fft)-1,abs(x_c1fft));
axis([0 31 0 20]);
title('三角波32点FFT')
for n=0:7
    if n<=3
        x_d(n+1)=4-n;
    end
    if n>=4
        x_d(n+1)=n-4;
    end
end
subplot(2,3,4)
stem(0:length(x_d)-1,x_d);
axis([0 8 0 5]);
title('反三角波')
subplot(2,3,5)
x_dfft=fft(x_d,8);
stem(0:length(x_dfft)-1,abs(x_dfft));
axis([0 8 0 20]);
title('反三角波8点FFT')
subplot(2,3,6)
x_d1fft=fft(x_d,32);
stem(0:length(x_d1fft)-1,abs(x_d1fft));
axis([0 31 0 20]);
title('反三角波32点FFT')
%% 保存图片
saveas(1, [Save_Path '三角波与反三角波'],'svg'); 
close all;
%% 第四小问
close all
N=16;
n=0:N-1;
df=1/16;
x=sin(2*pi*0.125*n)+cos(2*pi*(0.125+df)*n);
xfft=real(fft(x));
figure;
subplot(2,2,1)
stem(abs(xfft))
title('N=16,df=1/16')
subplot(2,2,2)
df=1/64;
x1=sin(2*pi*0.125*n)+cos(2*pi*(0.125+df)*n);
xfft1=real(fft(x1));
stem(abs(xfft1))
title('N=16,df=1/64')
N=128;
n=0:N-1;
df=1/16;
x2=sin(2*pi*0.125*n)+cos(2*pi*(0.125+df)*n);
xfft2=real(fft(x2));
subplot(2,2,3)
stem(abs(xfft2))
title('N=128,df=1/16')
df=1/64;
x3=sin(2*pi*0.125*n)+cos(2*pi*(0.125+df)*n);
xfft3=real(fft(x3));
subplot(2,2,4)
stem(abs(xfft3))
title('N=128,df=1/64')
%% 保存图片
saveas(1, [Save_Path '第四问'],'svg'); 
close all;
%% 第五问
close all
clc
%% FFT卷积
figure
n=0:15;
p=8;
q=2;
a=0.1;
f=0.0625;
x_a1=exp(-(n-p).^2./q);
x_b1=exp(-a*n).*sin(2*pi*f*n);
x_afft1=fft(x_a1,16);
x_bfft1=fft(x_b1,16);
c_conv=ifft(x_afft1.*x_bfft1);%循环卷积结果
subplot(2,1,1)
stem(0:length(c_conv)-1,c_conv)
axis([0 16 -1 2])
title('16点循环卷积')
subplot(2,1,2)
%c_linear_conv=conv(x_a1,x_b(1,:));
%利用fft计算线性卷积，不要用conv函数
x_afft1=fft(x_a1,32);
x_bfft1=fft(x_b1,32);
c_linear_conv=ifft(x_afft1.*x_bfft1);
stem(0:length(c_linear_conv)-1,c_linear_conv)
axis([0 31 -1 2])
title('线性卷积')
%% 保存图片
saveas(1, [Save_Path '第五问'],'svg'); 
close all;
%% 第六问
close all
clc
figure
rng(1000,"twister");%固定随机数
x_e=randn(1,512);%按行存储
x_efft=real(fft(x_e));
subplot(2,1,1)
stem(abs(x_efft))
title('卷积前频谱')
subplot(2,1,2)
x_econv=conv(x_e,x_c);
%一定要取实数！！！！！
x_econv_fft=real(fft(x_econv));
stem(abs(x_econv_fft))
title('线性卷积后频谱')
figure
subplot(2,1,1)
x_econvadd=overlap_add_convolution(x_e',x_c');
%一定要取实数！！！！！
x_econvadd_fft=real(fft(x_econvadd));
stem(abs(x_econvadd_fft))
title('重叠相加法卷积后频谱')
subplot(2,1,2)
x_econvkeep=overlap_keep_convolution(x_e',x_c');
%一定要取实数！！！！！
x_econvkeep_fft=real(fft(x_econvkeep));
stem(abs(x_econvkeep_fft))
title('重叠保留法卷积后频谱')
%% 保存图片
saveas(1, [Save_Path '第六问'],'svg'); 
saveas(2, [Save_Path '第六问2'],'svg'); 
close all;
%% 第七问
close all
clc
%% ===================【原理分析】===================
%% 知识点
%N点xn与yn补零到长度2N-1
%计算共轭的FFT
%Rxy=IDFT[X*(K)Y(K)]
%注意得到结果后需要重排顺序
%循环卷积与线性卷积的关系
%当L=N+M-1时两者相同，当L小，混叠，L大补零。
%%
close all
clc
figure
n=0:1:15;
X_a=exp(-(n-8).^2/2);        
X_b=exp(-(0.1).*n).*sin(2*pi*(0.0625)*n); 
k=length(X_a);
X_afft=fft(X_a,2*k);
X_bfft=fft(X_b,2*k);
rm=real(ifft(conj(X_afft).*X_bfft));
rm=[rm(k+2:2*k) rm(1:k)];
m=(-k+1):(k-1);
subplot(2,2,1)
stem(m,rm);
title('线性相关Xa,Xb')
%16点循环相关直接由线性相关延拓得到周期相关取主值区间得到
subplot(2,2,2)
rm1=[rm(16) rm(1:15)+rm(17:31)];
mm=0:15;
stem(mm,rm1);
title('循环相关Xa,Xb')
subplot(2,2,3)
rm=real(ifft(conj(X_bfft).*X_afft));
rm=[rm(k+2:2*k) rm(1:k)];
stem(m,rm);
title('线性相关Xb,Xa')
subplot(2,2,4)
rm1=[rm(16) rm(1:15)+rm(17:31)];
stem(mm,rm1);
title('循环相关Xb,Xa')
%% 保存图片
saveas(1, [Save_Path '第七问'],'svg'); 
close all;
%% 第八问
close all
clc
figure
X_afft=fft(X_a,2*k);
X_bfft=fft(X_b,2*k);
rm=real(ifft(conj(X_afft).*X_afft));
rm=[rm(k+2:2*k) rm(1:k)];
m=(-k+1):(k-1);
subplot(2,1,1)
stem(m,rm);
title('自相关Xa')
X_afft=fft(X_a,2*k);
X_bfft=fft(X_b,2*k);
rm=real(ifft(conj(X_bfft).*X_bfft));
rm=[rm(k+2:2*k) rm(1:k)];
m=(-k+1):(k-1);
subplot(2,1,2)
stem(m,rm);
title('自相关Xb')
%% 保存图片
saveas(1, [Save_Path '第八问'],'svg'); 
close all;








