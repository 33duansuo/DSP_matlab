clc
close all
clear all
%%
% FolderPath='E:\课程\海洋声信号实验\Data_12Week\Data_12Week';
% temp1=dir([FolderPath,'*.mat']);
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_CW_HFM.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_CW1.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_CW2.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_HFM.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_LFM.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\SndSig_CW.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\SndSig_CW_HFM.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\SndSig_HFM.mat"
%"E:\课程\海洋声信号实验\Data_12Week\Data_12Week\SndSig_LFM.mat"
Data=load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\SndSig_LFM.mat');
Data1=load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_LFM.mat');
% Data2=load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\RcvSig_CW2.mat');
Scv=Data.s_T;
True_Rcv=Data1.RcvSig;
% True_Rcv1=Data2.RcvSig;

fs=5000;
v=-3:0.5:3;
N=length(Scv);
tt = (0:N-1)/fs; % 时间向量
Rcv1=Signal_Receive(Scv,v(1),fs);
Rcv2=Signal_Receive(Scv,v(2),fs);
Rcv3=Signal_Receive(Scv,v(3),fs);
Rcv4=Signal_Receive(Scv,v(4),fs);
Rcv5=Signal_Receive(Scv,v(5),fs);
Rcv6=Signal_Receive(Scv,v(6),fs);
Rcv7=Signal_Receive(Scv,v(7),fs);
Rcv8=Signal_Receive(Scv,v(8),fs);
Rcv9=Signal_Receive(Scv,v(9),fs);
Rcv10=Signal_Receive(Scv,v(10),fs);
Rcv11=Signal_Receive(Scv,v(11),fs);
Rcv12=Signal_Receive(Scv,v(12),fs);
Rcv13=Signal_Receive(Scv,v(13),fs);

% %% 匹配滤波1
% % 匹配滤波器（脉冲响应为信号的时间反转）
% matched_filter = flipud(Scv); 
% figure
% % 进行匹配滤波
% Cwoutput_signal = conv(matched_filter,Rcv1,'same');
% plot(tt,Cwoutput_signal);

%% 匹配滤波2
%相关函数卷积结果
matched_filter = flipud(True_Rcv);
% matched_filter1 = flipud(True_Rcv1);
% 使用 xcorr 计算自相关
%[Cwoutput, lags] = xcorr(Input_signal,Input_signal);%如果s与noise不相关可以直接用
for i=1:1:13
    subplot(4,4,i)
    Cwoutput=conv(eval(['Rcv' num2str(i)]), matched_filter);
    %[Cwoutput, lags] = xcorr(True_Rcv,eval(['Rcv' num2str(i)]));
    % 为了绘制输出信号，生成相应的时间轴
    lag_time = (0:length(Cwoutput)-1)/fs;
    plot(lag_time,Cwoutput)
    [peek(i),~]=max(Cwoutput);
    if i==13
        [~,max_time1]=max(Cwoutput);
        Aim_Length=1500*(max_time1/fs-N/fs)/2;
        title(['LFM目标距离为' num2str(Aim_Length) 'm'])
    end
end
[~,location]=max(peek);
V_Estimate=v(location);

% figure
% for i=1:1:13
%     subplot(4,4,i)
%     Cwoutput=conv(eval(['Rcv' num2str(i)]), matched_filter1);
%     %[Cwoutput, lags1] = xcorr(True_Rcv1,eval(['Rcv' num2str(i)]));
%     % 为了绘制输出信号，生成相应的时间轴
%     lag_time = (0:length(Cwoutput)-1)/fs;
%     plot(lag_time,Cwoutput)
%     [peek1(i),~]=max(Cwoutput);
%     if i==12
%         [~,max_time]=max(Cwoutput);
%         Aim_Length1=1500*(max_time/fs-N/fs)/2;
%         title(['目标距离为' num2str(Aim_Length1) 'm'])
%     end
% end
% [~,location]=max(peek1);
% V_Estimate1=v(location);

disp(['====== LFM信号速度估计 ',num2str(V_Estimate),'m/s ======'])
disp(['====== LFM信号距离估计 ',num2str(Aim_Length),'m ======'])
% disp(['====== CW2信号速度估计 ',num2str(V_Estimate1),'m/s ======'])
% disp(['====== CW2信号距离估计 ',num2str(Aim_Length1),'m ======'])

%%
clc
close all
clear all
Data=load('E:\课程\海洋声信号实验\Data_12Week\Data_12Week\arrayData.mat');
xxdelay=Data.arrayData;
%% 阵列处理

close all;
M=16;
Fs=5000;
c=1500;
% d=c/Fs*0.5;
d=2;
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
Ypp_Energy=normalize(Ypp_Energy,1,"norm",5);
imagesc(theta1*180/pi,tt,(Ypp_Energy'));
%%
figure(3)
plot(theta1/pi*180,10*log10(BeamEnergy/max(BeamEnergy)),'b-','LineWidth',1)
xlabel('\fontsize{15}angle (°)')
ylabel('\fontsize{15}Power (dB)')
title('波束扫描结果')
[~,Aim_angle]=max(BeamEnergy);
disp(['====== 目标方位估计为 ',num2str(Aim_angle-1),'° ======'])






%%
function  s_R=Signal_Receive(s_T,v,fs)
    T=length(s_T)/fs;
    N=T*fs;
    filterOrder=10;
    c=1500;
    delta=2*v/c;
    RcvSigLen=round((1-delta)*N);
    s_R=zeros(1,RcvSigLen);
    ii=1;
    dFS=1/(1-delta)*ii;
    DI=round(dFS);%数字时延中的整数部分
    while(DI<=N)
        taof=-(dFS-DI);%数字时延中的小数部分
        if DI < filterOrder/2+1
            s_R(ii)=newDelayFilter([zeros(1,filterOrder/2+1-DI) s_T(1:DI) s_T(ii+DI+(1:filterOrder/2))], -taof, filterOrder/2); 
        elseif DI > N-filterOrder/2
            s_R(ii)=newDelayFilter([s_T(DI+(-filterOrder/2:-1)) s_T(DI:N) zeros(1,filterOrder/2+1-(N-DI+1))], -taof, filterOrder/2); 
        else
            s_R(ii)=newDelayFilter([s_T(DI+(-filterOrder/2:filterOrder/2))], -taof, filterOrder/2); 
        end
        ii=ii+1;
        dFS=1/(1-delta)*ii;
        DI=round(dFS);%数字时延中的整数部分
    end
    % rcVnn=0:RcvSigLen-1;
    % rcVtt=rcVnn/fs;
    
    % figure(1)
    % plot(tt,s_T,'b-')
    % hold on
    % plot(rcVtt,s_R,'r-')
    % xlabel('\fontsize{15} 时间/s')
    % ylabel('\fontsize{15} 幅度')
    % grid on
    % legend('\fontsize{15}发射信号',' \fontsize{15}接收信号')
    % title('\fontsize{15}收发信号时域波形')
    
    % 
    % TransmitDataFFT=fft(s_T,N);
    % TransmitDataPsd=((abs(TransmitDataFFT)).^2)/N/fs;
    % 
    % deltaf=fs/N;
    % ff=0:deltaf:fs-deltaf;
    % 
    % if RcvSigLen>N
    %     RcvData=s_R(1:N);
    % else
    %     RcvData=[s_R zeros(1,N-RcvSigLen)];
    % end
    % RcvDataFFT=fft(RcvData,N);
    % RcvDataPsd=((abs(RcvDataFFT)).^2)/N/fs;
    
    % figure(2)
    % plot(ff,10*log10(TransmitDataPsd),'b-')
    % hold on
    % plot(ff,10*log10(RcvDataPsd),'r-')
    % 
    % legend('\fontsize{15}发射信号功率谱','\fontsize{15}接收信号功率谱')
    % xlabel('\fontsize{15}频率/Hz')
    % ylabel('\fontsize{15}功率谱/dB')
    % grid on
    % axis([f0-10 f0+30 -350 0])
    % title('\fontsize{15}收发信号功率谱')
end