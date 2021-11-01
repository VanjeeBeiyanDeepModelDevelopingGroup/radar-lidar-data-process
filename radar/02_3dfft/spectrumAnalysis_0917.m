function [rangeSpectrum,angleSpectrum] = spectrumAnalysis_0917(dataCube3dFFT)
[a,b,c] = size(dataCube3dFFT);
%spectrumAnalysis 频谱分析，画出一维波形，包括距离的和角度的
% rangeDoppler_sum = squeeze(sum(abs(dataCube3dFFT),2));
% rangeSpectrum = 20*log10(squeeze(sum(rangeDoppler_sum,2)));
rangeDoppler_sum = squeeze(max(abs(dataCube3dFFT),[],2));
rangeSpectrum = 20*log10(squeeze(max(rangeDoppler_sum,[],2)));
figure(4);
subplot(4,2,1);plot(rangeSpectrum);
% 在range 图上CFAR提取再做角度频谱
angleSpectrum=[];
%% 短时傅里叶变换
range_tempArray = squeeze(dataCube3dFFT(:,1,1));
subplot(4,2,2);
stft(range_tempArray,a,'Window',kaiser(a,5),'OverlapLength',4,'FFTLength',a);title('距离维度的STFT');
vel_tempArray = squeeze(dataCube3dFFT(1,1,:));
subplot(4,2,3);
stft(vel_tempArray,c,'Window',kaiser(c,5),'OverlapLength',2,'FFTLength',c);title('速度维度的STFT');

end

