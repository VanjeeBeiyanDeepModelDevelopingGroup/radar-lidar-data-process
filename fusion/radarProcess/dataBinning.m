function [rangeDoppler_sum,rangeAngle_sum,DopplerAngle_sum] = dataBinning(dataCube3dFFT)
%%dataBinning 数据累积，谱分析等方法，为了提升信噪比，压旁瓣等

[range_usrr,angle_usrr,doppler_usrr] = size(dataCube3dFFT);
% rangeDoppler_sum = sum(abs(dataCube3dFFT),2);
%% maxpool：取最大值
% rangeDoppler_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],2)));
% rangeAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],3)));
% DopplerAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],1)));
%% 非相干累积：幅值累加
% rangeDoppler_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),2)));
% rangeAngle_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),3)));
% DopplerAngle_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),1)));
%% 相干累积：复数相加
% rangeDoppler_sum = 20*log10(squeeze(abs(sum(dataCube3dFFT,2))));
% rangeAngle_sum = 20*log10(squeeze(abs(sum(dataCube3dFFT,3))));
% DopplerAngle_sum = 20*log10(squeeze(abs(sum(dataCube3dFFT,1))));
%% 非相干累乘：幅值累乘
% abs_dataCube3dFFT = abs(dataCube3dFFT);
% % 这里需要归一化，否则就太大了
% max_datacube = max(abs_dataCube3dFFT(:));
% min_datacube = min(abs_dataCube3dFFT(:));
% norm_abs_dataCube3dFFT = (abs_dataCube3dFFT-min_datacube)/(max_datacube-min_datacube)+1.0;
% rangeDoppler_sum = ones([range_usrr,doppler_usrr]);
% for i=1:angle_usrr
%     rangeDoppler_sum = rangeDoppler_sum.*squeeze(norm_abs_dataCube3dFFT(:,i,:));
% end
% % rangeDoppler_sum = 20*log10(rangeDoppler_sum);
% rangeAngle_sum = ones([range_usrr,angle_usrr]);
% for i=1:doppler_usrr
%     rangeAngle_sum = rangeAngle_sum.*(norm_abs_dataCube3dFFT(:,:,i));
% end
% 
% DopplerAngle_sum = ones(size(doppler_usrr,angle_usrr));
% for i=1:range_usrr
%     DopplerAngle_sum = DopplerAngle_sum.*squeeze(norm_abs_dataCube3dFFT(i,:,:));
% end
%% RODNET方法：取有相对速度的进行累加
% RODNET方法仅针对Range-Angle Map，与其他的热力图无关
% 这里其他两个维度的热力图暂时采用取最大值的方式
% 第1维 range   512
% 第2维 angle   64（只有12,补零到64）
% 第3维 doppler 32（32个chirp）
% 这个维数是怎么回事？？？
% 先把速度为0的部分筛掉
rangeDoppler_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],2)));
%___
R_D = squeeze(max(abs(dataCube3dFFT),[],2));
R_D_2 = 20*log10(R_D);
% figure(6);
% subplot(2,1,1);
% mesh(R_D);
% xlabel('doppler');
% ylabel('range');
% subplot(2,1,2);
% mesh(R_D_2);
% xlabel('doppler');
% ylabel('range');
%___
rangeAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],3)));
%___
R_A = squeeze(max(abs(dataCube3dFFT),[],3));
R_A_2 = 20*log10(R_A);
% figure(7);
% subplot(2,1,1);
% mesh(R_A);
% xlabel('angle');
% ylabel('range');
% subplot(2,1,2);
% mesh(R_A_2);
% xlabel('angle');
% ylabel('range');
%___
DopplerAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],1)));

%% 平台直方图对比度增强CLAHE
% enhancedImg = adapthisteq(rangeDoppler_sum);% ,'clipLimit',0.1,'Distribution','rayleigh');
% figure;subplot(1,2,1);imagesc(rangeDoppler_sum);
% subplot(1,2,2);imagesc(enhancedImg);
end

