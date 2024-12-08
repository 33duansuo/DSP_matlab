clc
close all
clear
Fs=5000;
v=5;
c=1500;
Tao=2;
t=0:1/Fs:Tao-1/Fs;
N=2000;
tx=chirp(t,300,2,300);
Scale_change_factor=2*v/c;
zerolength=ceil(Scale_change_factor*Tao*Fs);
rx=zeros(1,2000+zerolength*2);
for i=1:2000
    DFS=t(i)*(Scale_change_factor)*Fs;
    BUFF=Delayfilter([zeros(1,i-1) tx(i) zeros(1,2000-i)],DFS,11,zerolength);
    rx=rx+BUFF;
end
ffts=fft(tx,Fs);
ff=(0:Fs-1);
plot(ff,10*log10(abs(ffts.^2)/N/Fs)+120);
axis([0 2500 0 120])
figure
ffts1=fft(rx,Fs);
plot(ff,10*log10(abs(ffts1.^2)/N/Fs)+120);
axis([0 2500 0 120])
% spectrogram(tx, hamming(1024), 120, 1024, Fs, 'yaxis');
% figure
% spectrogram(rx, hamming(1024), 120, 1024, Fs, 'yaxis');
