function outdata=Delayfilter(inputdata,DFS,N,zerolength)
%(inputdata,DFS,N,zerolength)
%DFS为数字延时
%N为需要的滤波器阶数
%zerolength为输出数据补零长度

%先整体延时
DI=round(DFS);%数字整数时延
tao=DFS-DI;%数字小数时延应该在-0.5到0.5之间
num1=zerolength+DI;
num2=zerolength-DI;
   xx=[zeros(1,num1) inputdata zeros(1,num2)];
M=N/2;
if tao~=0
    k=-M:M;
    hd=sinc(k-tao);
    hnWin=hanning(2*M+1)';
    hd=hd.*hnWin;
    h=hd/sum(hd);
    xxconv=conv(xx,h);
    outdata=xxconv(2*M+1:length(xx)+2*M);%取完全滤波数据
end
if tao==0
   outdata=xx;
end


