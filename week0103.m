clc
clear
close all
%%
Data = load('E:\课程\海洋声信号实验\AnechoicTankExperimentData\AnechoicTankExperimentData\PassiveArrayRcvData_ShipNoise_1_2.mat');
%%
TraceBeamData = Data.ArrayData;
%% 阵列处理

close all;
M=6;
Fs=5000;
c=1500;
d=0.09;
filterOrder=10;
zeroslength=ceil((M-1)*d)/c*Fs;
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
        Yp=Yp+outData;
    end
    Yp=Yp/M;
    Ypp(ii,:)=Yp;
    YpFFT(ii,:)=fft(Yp)/WL;
    BeamEnergy(ii)=sum(abs(YpFFT(ii,:)).^2,'all');
 end
%%
close all
figure(2)
tt=(0:length(YpFFT)-1)/Fs;
for i=1:181
    Ypp_Energy(i,:)=cumsum(abs(Ypp(i,:)).^2);%能量累加
end


%单帧为1s
for i=1:181
    for tt_1s=2:1:200
        Ypp_Energy_1s(i,tt_1s)=Ypp_Energy(i,tt_1s*5000)-Ypp_Energy(i,(tt_1s-1)*5000);
    end
end

%7帧叠加做估计
for tt_1s=2:1:200
    if tt>=4&tt<=197
        Ypp_Energy_1s(i,tt_1s)=sum(Ypp_Energy_1s(i,tt_1s-3:tt_1s+3));
    end
    if tt<4
        Ypp_Energy_1s(i,tt_1s)=sum(Ypp_Energy_1s(i,1:tt_1s));
    end
    if tt>197
        Ypp_Energy_1s(i,tt_1s)=sum(Ypp_Energy_1s(i,tt_1s:end));
    end
end

Ypp_Energy_1s=normalize(Ypp_Energy_1s,1,"norm",5);
Ypp_Energy_1s(:,1)=[];%删除第一行

%% 作图
close all
tt_1s=2:1:200;
imagesc(theta1*180/pi,tt_1s,(Ypp_Energy_1s'));
colormap jet;  % 设置颜色映射
colorbar;  % 显示颜色条
[~,Aim_angle10s]=max(Ypp_Energy_1s(:,9));
[~,Aim_angle50s]=max(Ypp_Energy_1s(:,49));
[~,Aim_angle100s]=max(Ypp_Energy_1s(:,99));
[~,Aim_angle150s]=max(Ypp_Energy_1s(:,149));
[~,Aim_angle200s]=max(Ypp_Energy_1s(:,199));
text(theta1(Aim_angle10s)/pi*180,9,['10s时刻方位' num2str(Aim_angle10s) '°'])%\leftarrow表示左箭头
text(theta1(Aim_angle50s)/pi*180,49,['50s时刻方位' num2str(Aim_angle50s) '°'])%\leftarrow表示左箭头
text(theta1(Aim_angle100s)/pi*180,99,['100s时刻方位' num2str(Aim_angle100s) '°'])%\leftarrow表示左箭头
text(theta1(Aim_angle150s)/pi*180,149,['150s时刻方位' num2str(Aim_angle150s) '°'])%\leftarrow表示左箭头
text(theta1(Aim_angle200s)/pi*180,199,['200s时刻方位' num2str(Aim_angle200s) '°'])%\leftarrow表示左箭头


%%
figure(3)
plot(theta1/pi*180,10*log10(BeamEnergy/max(BeamEnergy)),'b-','LineWidth',1)
xlabel('\fontsize{15}angle (°)')
ylabel('\fontsize{15}Power (dB)')
title('波束扫描结果')
[~,Aim_angle]=max(BeamEnergy);
disp(['====== 目标方位估计为 ',num2str(Aim_angle-1),'° ======'])