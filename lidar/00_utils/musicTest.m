clear
clc
% load('datacube.mat');
% [range,doppler,angle] = size(datacube);
load('Temp1_');
%% 验证多次fft跟fft2的区别：没有区别
% for i=1:angle
%     rd_plane = squeeze(datacube(:,:,i));
%     fft_rd_plane = fft2(rd_plane);
%     figure;imagesc(abs(fft_rd_plane));
% end

% fft_rd_map = fft(datacube,range,1);
% fft_rd_map = fft(fft_rd_map,doppler,2);
% fft_3dcube = fft(fft_rd_map,angle,3);
% rd_map = sum(abs(fft_3dcube(:,:,1:8)),3);
% figure;imagesc(rd_map);

% x = datacube(:,2,2);
x=Temp1;
fs = length(x);
nwin = fs;
num=floor(fs/2)-1;
figure(4);clf();subplot(4,1,1);plot(imag(x));hold on;
plot(real(x));legend('虚部','实部');title('时域信号');
% %% 给信号加窗
% % 加窗
angleBin_num=64;
angleWin = hanning(angleBin_num);  
angleWin = angleWin(1: (angleBin_num / 2));
angleWinLen               = length(angleWin);
angleWindowCoeffVec       = ones(angleBin_num, 1);
angleWindowCoeffVec(1:angleWinLen) = angleWin;
angleWindowCoeffVec(angleBin_num-angleWinLen+1:angleBin_num) = angleWindowCoeffVec(angleWinLen:-1:1);
angleWin = angleWindowCoeffVec;
Temp1(1:angleBin_num) = bsxfun(@times, Temp1(1:angleBin_num), angleWin.');
x2=Temp1;
% subplot(4,1,2);plot(imag(x2));hold on;
% plot(real(x2));legend('虚部','实部');title('时域信号');
% subplot(4,1,3);plot(angleWin);
%% fft
fft_x = fft(x);
fft_x2 = fft(x2);
subplot(4,1,2);plot(20*log10(abs(fft_x)));title('fft');hold on
plot(20*log10(abs(fft_x2)));hold off;legend('加窗前','加窗后');
%% music方法
for i=2:num
[P,f] = pmusic(x,i,fs,fs); % input,p,nfft,fs,nwin,overlap
[P2,f2] = pmusic(x,[Inf,1.1],fs,fs,i); % Window length = 7
subplot(4,1,3);plot(f,20*log10(abs(P)));
str = ['MUSIC-功率谱随子空间维度增大，空间维度：',num2str(i)];
title(str);grid on
subplot(4,1,4);plot(f2,20*log10(abs(P2)))
title('MUSIC-功率谱随nwin增大');
xlabel('freq');ylabel('dim');
xlabel 'Frequency (Hz)'

pause(0.01);
end
