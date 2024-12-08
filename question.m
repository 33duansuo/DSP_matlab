clc;
clear all;
close all;
Fs = 8000;         
T = 2;           
N = Fs * T;      
tt = (0:N-1) / Fs;
F1=350;
F2=750;
s =randn(N, 1);
%% 带通滤波器设计
fl_fir=375;
fh_fir=750;
a_fir=1;
b_fir=fir1(64,[fl_fir/(Fs/2),fh_fir/(Fs/2)]);
%% 阵列延时
close all;
d=1;
theta0=30.98/180*pi;
c=1500;
MM=16;
addlength=80;
filterOrder=10;
for m=1:MM
tm=(m-1)*d*cos(theta0)/c*Fs;
outdata=Delayfilter(s',tm,filterOrder,addlength);
arraydelay(m,:)=outdata;
end

%% 波束形成
close all;
M=16;
d=1;
c=1500;
wave_length=1;%波束扫描间隔
zeroslength=(M-1)*d/c*Fs;
theta1=(0:wave_length:180)*pi/180;
BeamNum=length(theta1);
BeamEnergy=zeros(1,BeamNum);
WL=length(arraydelay)+2*zeroslength;
 for ii=1:BeamNum%波束
     Yp=zeros(1,WL);%预成波束数据
    for jj=1:M%阵元
        DFS=-(jj-1)*d*cos(theta1(ii))/c*Fs;
        DI=DFS-round(DFS);
        outData=Delayfilter(arraydelay(jj,:),DFS,filterOrder,zeroslength);
        % outData=filter(b_fir,a_fir,outData);
        Yp=Yp+outData;
    end
    Yp=Yp/M;
     YpFFT=fft(Yp)/WL;
    BeamEnergy(ii)=sum(abs(YpFFT).^2,'all');
 end  
figure(2)
plot(theta1/pi*180,10*log10(BeamEnergy/max(BeamEnergy)),'b-','LineWidth',1)
xlabel('\fontsize{15}angle (°)')
ylabel('\fontsize{15}Power (dB)')

%% 方向估计
[~,locs]=findpeaks(BeamEnergy,'NPeaks',1,'SortStr','descend');
if length(locs)==1
    maxValueIdx=locs;
else
   [~,maxValueIdx]=max(BeamEnergy);
end
Z0=sqrt(BeamEnergy(maxValueIdx-1));
Z1=sqrt(BeamEnergy(maxValueIdx));
Z2=sqrt(BeamEnergy(maxValueIdx+1));
modifyValue=(Z2-Z0)/(4*Z1-2*Z0-2*Z2);
if abs(modifyValue)>1
    modifyValue=sign(modifyValue);
end
estang=(maxValueIdx-1)*wave_length+modifyValue*wave_length;
disp(['====== 真实目标方位',num2str(theta0*180/pi),'° ======'])
disp(['====== 方向估计完成  目标方位',num2str(estang),'° ======'])

















