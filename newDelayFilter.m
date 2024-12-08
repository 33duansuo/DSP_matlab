%%��ȷС��ʱ���˲���
% inData:����ο�����   ����ΪFrameLen+2*M
% taof:С������ʱ����
% M:2*MΪС��ʱ���˲����Ľ���
% outData:���ʱ�Ӻ�����У�һ��
%�˲�������ΪN=2*M+1��Ⱥʱ��Ϊ(N-1)/2=M  
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