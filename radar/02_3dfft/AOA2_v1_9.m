% 取range下的slow_time*channel矩阵输入music算法计算angle谱响应
% 特点：遍历按range列排序后的cfar结果，构造channel*slow time的matrix进行music谱估计
function [ spectrumAngle, pointRAarr ] = AOA2_v1_9( dopplerOut, cfarOut, TX_num, RX_num, numADCSamples, dopplerBin_num)
    if isempty(cfarOut)
        spectrumAngle = [];
        pointRAarr = zeros(numADCSamples,1);
        return
    end
    % angleBin_num 输入真实天线个数
    L = 360;
    if TX_num == 3
        angleNum = (TX_num - 1) * RX_num;
    else
        angleNum = TX_num * RX_num;
    end
    dopplerOut3d = zeros(numADCSamples, dopplerBin_num,angleNum);
    spectrumAngle = zeros(numADCSamples,L);
    % 保存目标点
    pointRAarr = zeros(numADCSamples,1);
    signalMat = [];
    lastRangeIndx = cfarOut(1,1); 
    for r = 1: length(cfarOut(:,1))
         rangeIndx = cfarOut(r,1);
         dopplerIndx = cfarOut(r,2);
         channelSignal = reshape(dopplerOut(rangeIndx,dopplerIndx,:,:),[1,TX_num*RX_num]);
         channelSignal = channelSignal(1:angleNum);
         if rangeIndx == lastRangeIndx
             signalMat = [signalMat,channelSignal'];
         else
             temp = musicAlg(signalMat,L);
             spectrumAngle(rangeIndx,:) = temp;
             % 角度最大值
             [angleMax, angleIndx] = max(temp);
             pointRAarr(rangeIndx,1) = angleIndx;
             signalMat = [];
             signalMat = [signalMat,channelSignal'];
         end
         lastRangeIndx = rangeIndx;
    end
end