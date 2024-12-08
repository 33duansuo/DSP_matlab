%%精确小数时延滤波器
% inData:输入参考序列   长度为FrameLen+2*M
% taof:小数数字时延量
% M:2*M为小数时延滤波器的阶数
% outData:输出时延后的序列，一行
%滤波器长度为N=2*M+1，群时延为(N-1)/2=M  
function outData =newDelayFilter(inData, taof, M) 
    [sa,~]=size(inData);
    if sa>1
        xx=inData.';
    else
        xx=inData;
    end
   
    if abs(taof)<0.000001
        yy1=xx;
        yy=yy1(M+1);
    else
        k=-M:M;
        hd=sinc(k-taof);
        hnWin=hanning(2*M+1)';
        hd=hd.*hnWin;
        h=hd/sum(hd);
        yy=xx*h';
    end
    
    if sa>1
        outData=yy';
    else
        outData=yy;
    end

end