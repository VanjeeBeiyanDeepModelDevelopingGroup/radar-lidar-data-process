function [rangeFFTwaveCluster] = splitRangeFFT(rangeFFT,cfarResult)
%   �˺�����RangeFFT����CFAR�Ľ����ֿ���ÿ����ֵ����Ӧ��һ��rangeFFT����
%   �ò��ν�����theta����ϣ��������ά�ռ��ϵ�rangeFFT����
%   ��ְ�������CFAR��index�������ֵ��������ķ�ֵ��⣬��ȡCFAR�����ķǷ�ֵ��Ϊ�������ֵ���ֵ��
% indexList =
figure;plot(rangeFFT);hold on
cnt = 0;
for i=1:length(cfarResult)
    cnt = cnt+1;
    if(cfarResult(i)==1)
        fprintf("%d\n", cnt);
        plot(i,rangeFFT(i),'ro'); hold on;
    end
end

% figure;plot(cfarResult); hold on


cfarResult = sort(cfarResult);
cfarLength = length(cfarResult);
[rangeFFTrows,rangeFFTcols] = size(rangeFFT);
rangeFFTwaveCluster = zeros(cfarLength,rangeFFTcols);
startIndex = 0;
for i=1:cfarLength-1
    cfarResultIndex = cfarResult(i);
    cfarResultIndexP1 = cfarResult(i+1);
    cutIndex = floor(cfarResultIndex+cfarResultIndexP1)/2;
    rangeFFTwaveCluster(i,startIndex:cutIndex) = rangeFFT(startIndex:cutIndex);
    startIndex = cutIndex;
end
% ��ʾ��άͼ��
[X,Y]=meshgrid(rangeFFTwaveCluster);
figure;surf(X,Y,rangeFFTwaveCluster);
