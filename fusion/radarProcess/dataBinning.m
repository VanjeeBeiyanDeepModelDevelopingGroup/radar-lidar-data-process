function [rangeDoppler_sum,rangeAngle_sum,DopplerAngle_sum] = dataBinning(dataCube3dFFT)
%%dataBinning 数据累积，谱分析等方法，为了提升信噪比，压旁瓣等

[range_usrr,angle_usrr,doppler_usrr] = size(dataCube3dFFT);
x=1:range_usrr;
y=1:angle_usrr;
[X,Y]=meshgrid(x,y);

% rangeDoppler_sum = sum(abs(dataCube3dFFT),2);
%% 查看所有速度页下的range angle图
% figure(4);
% rows = 4;
% cols = doppler_usrr/rows;
% for i =1:doppler_usrr
%     subplot(rows,cols,i);% imagesc(abs(dataCube3dFFT(:,:,i)));
%     % 画3d surf图
%     figure;surf(X,Y,abs(dataCube3dFFT(:,:,i-4))');
%     set(gca,'YDIR','normal');
% end
%% maxpool：取最大值
% rangeDoppler_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],2)));
% rangeAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],3)));
% DopplerAngle_sum = 20*log10(squeeze(max(abs(dataCube3dFFT),[],1)));
%% 非相干累积：幅值累加
% rangeDoppler_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),2)));
% rangeAngle_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),3)));
% DopplerAngle_sum = 20*log10(squeeze(sum(abs(dataCube3dFFT),1)));

rangeDoppler_sum = (squeeze(sum(abs(dataCube3dFFT),2)));
rangeAngle_sum = (squeeze(sum(abs(dataCube3dFFT),3)));
DopplerAngle_sum = (squeeze(sum(abs(dataCube3dFFT),1)));
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


%% 平台直方图对比度增强CLAHE
% enhancedImg = adapthisteq(rangeDoppler_sum);% ,'clipLimit',0.1,'Distribution','rayleigh');
% figure;subplot(1,2,1);imagesc(rangeDoppler_sum);
% subplot(1,2,2);imagesc(enhancedImg);
%% 归一化
rangeDoppler_sum = (rangeDoppler_sum-min(rangeDoppler_sum(:)))/(max(rangeDoppler_sum(:))-min(rangeDoppler_sum(:)));
rangeAngle_sum = (rangeAngle_sum-min(rangeAngle_sum(:)))/(max(rangeAngle_sum(:))-min(rangeAngle_sum(:)));
DopplerAngle_sum = (DopplerAngle_sum-min(DopplerAngle_sum(:)))/(max(DopplerAngle_sum(:))-min(DopplerAngle_sum(:)));
end

