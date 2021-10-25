function [rangeFFTwaveCluster] = splitRangeFFT(rangeFFT,cfarResult)
%   此函数将RangeFFT按照CFAR的结果拆分开，每个峰值都对应了一个rangeFFT波形
%   该波形将与其theta角组合，组成在三维空间上的rangeFFT波形
%   拆分按照输入CFAR的index，将最大值和其邻域的幅值拆解，新取CFAR滑窗的非峰值作为其他部分的数值。
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
% 显示二维图像
[X,Y]=meshgrid(rangeFFTwaveCluster);
figure;surf(X,Y,rangeFFTwaveCluster);
