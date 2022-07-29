%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults function                                              %%%
%%%     radarDataAll- 毫米波雷达数据                                        %%%
%%%     lidarDataFrame- 激光雷达数据                                        %%%
%%%     frameIndex- 想查看的帧序列位置                                      %%%
%%%     methodSign- 具体采用的增加最大速度的标识位                          %%%
%%%                                                                         %%%
%%%     distanceCoor- 横坐标距离刻度                                        %%%
%%%     velocityCoor- 纵坐标速度刻度                                        %%%
%%%     distanceOut- 选定帧的CFAR后筛的点的距离信息                         %%%
%%%     velocityOut- 选定帧的CFAR后筛的点的速度信息                         %%%
%%%     FinalResult- 选定帧的CFAR后筛的点的距离，速度，xyz轴坐标信息        %%%
%%%     mmwavedata- 一维FFT结果的绝对值                                     %%%
%%%     dopplerSum- doppler fft结果的加和值                                 %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.05.12 version 1.3                            %%%
%%%     修改：增加多种方法对速度的处理                                      %%%
%%%           方法选择用methodSign标志位进行区别                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [ distanceCoor, velocityCoor, distanceOut, velocityOut, FinalResult, mmwavedata,dopplerSum,pcStrc ] = mmwaveResults_v1_3(radarDataAll, lidarDataFrame, frameIndex, methodSign, lineId)
function [ distanceCoor, velocityCoor, distanceOut, velocityOut, FinalResult, mmwavedata,dopplerSum,pcStrc ] = mmwaveResults_v1_3(radarDataAll, frameIndex, methodSign, lineId)

numADCSamples = param.numADCSamples;
% numChirpsPerFrame = param.numChirpsPerFrame;
RX_num = param.RX_num;
TX_num = param.TX_num;
% numDataPerFrame = numChirpsPerFrame * RX_num * numADCSamples * 2;

frame_index = str2num(frameIndex);

% data_num = 128
[data_num,~,~,~] = size(radarDataAll);
dopplerBin_num = data_num/RX_num;
% data_num

%%% range FFT %%%
rangeWin = hanning(numADCSamples);
rangeWin = rangeWin(1: (numADCSamples / 2));
RangeWinLen               = length(rangeWin);
RangeWindowCoeffVec       = ones(numADCSamples, 1);
RangeWindowCoeffVec(1:RangeWinLen) = rangeWin;
RangeWindowCoeffVec(numADCSamples-RangeWinLen+1:numADCSamples) = RangeWindowCoeffVec(RangeWinLen:-1:1);
rangeWin   = RangeWindowCoeffVec;
% rangeWin_size = size(rangeWin.')
% 汉宁窗是512维
rangeOut = bsxfun(@times, radarDataAll, rangeWin.');
% rangeOut_size = size(rangeOut)
% 所有数据构成128x512x3x8，即32x4个chirp，512个采样点，3个发送天线，8帧
% 沿着第二个维度，做一维fft
rangeOut = fft(rangeOut, numADCSamples, 2);
% figure;
% rangeOutNew = zeros(size(rangeOut));
% [chirps,ads,txs,frames] = size(rangeOut);
% for chirp=1:chirps
%     
% [~,rangeOut] = pmusic(rangeOut,64,512,Fs,'whole');

c = param.c;
digOutSampleRate = param.digOutSampleRate;
% freqSlopeConst = 56.25;
freqSlopeConst = param.freqSlopeConst;  % 测试那个指令
slope = freqSlopeConst * 1e12;

rangeAbs = abs(rangeOut);

distanceCoor = ((1: numADCSamples) * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
% 返回的是第n帧的第1个天线的第一个chirp的距离FFT
% 更改，返回当前帧的第一个天线的所有chirp的range-FFT
mmwavedata = rangeAbs(:, :, 1, frame_index);
% mmwavedataSize = size(mmwavedata)
% 尝试一下拆分
% guardLen = 2;
% winLen = 16;
% thresholdScale = 40;
% noiseDivShift = ceil(log2(2 * winLen));
% cfarResult = cfar_ca( mmwavedata(1,:), numADCSamples, thresholdScale, noiseDivShift, guardLen, winLen);
% rangeFFTwaveCluster = splitRangeFFT(mmwavedata(1,:),cfarResult);

%%% doppler FFT %%%
dopplerWin = hanning(dopplerBin_num);  % 32个chirp
dopplerWin = dopplerWin(1: (dopplerBin_num / 2));
dopplerWinLen               = length(dopplerWin);
dopplerWindowCoeffVec       = ones(dopplerBin_num, 1);
dopplerWindowCoeffVec(1:dopplerWinLen) = dopplerWin;
dopplerWindowCoeffVec(dopplerBin_num-dopplerWinLen+1:dopplerBin_num) = dopplerWindowCoeffVec(dopplerWinLen:-1:1);
dopplerWin   = dopplerWindowCoeffVec;

dopplerIn = zeros(numADCSamples, dopplerBin_num, RX_num, TX_num);% 512x32x4x3

for n = 1: numADCSamples

    for i = 1: TX_num

        for m = 1: RX_num

            dopplerIn(n, :, m, i) = rangeOut((0: (dopplerBin_num - 1)) * RX_num + m, n, i, frame_index);

        end

    end

end
dopplerIn = bsxfun(@times, dopplerIn, dopplerWin.');
% dopplerInSize = size(dopplerIn)
dopplerOut = fftshift(fft(dopplerIn, dopplerBin_num, 2), 2);
% save dopplerOut dopplerOut
% dopplerOutSize = size(dopplerOut)
dopplerLog2Abs = log2(abs(dopplerOut));
% 画出进行angle-fft之前的xcube
% figure(2);
% dopplerLog2Abs_3d(:,:,1:4) = dopplerLog2Abs(:,:,:,1);
% dopplerLog2Abs_3d(:,:,5:8) = dopplerLog2Abs(:,:,:,2);
% dopplerLog2Abs_3d(:,:,9:12) = dopplerLog2Abs(:,:,:,3);
% for i=1:12
%     subplot(3,4,i);
%     imagesc(dopplerLog2Abs_3d(:,:,i));set(gca,'YDir','normal');title(num2str(i));
% end
% dopplerLog2AbsSize = size(dopplerLog2Abs)
dopplerSum = sum(dopplerLog2Abs, [3 4]);% 这里应该做beamforming吧，可以做beamforming，可以试试
% dopplerSumSize = size(dopplerSum)
dopplerSum = squeeze(dopplerSum); 
% dopplerSize = size(dopplerSum)

startFreqConst = param.startFreqConst;
startFreq = startFreqConst * 1e9;
bandwidth = (slope * numADCSamples) / (digOutSampleRate * 1e3);
adcStartTimeConst = param.adcStartTimeConst;
adcStart = adcStartTimeConst * 1e-6;
centerFreq = startFreq + bandwidth * 0.5 + adcStart * slope;
idleTimeConst = param.idleTimeConst;
rampEndTime = param.rampEndTime;
chirpInterval = (idleTimeConst + rampEndTime) * 1e-6;

velocityCoor = (((- dopplerBin_num / 2): (dopplerBin_num / 2 - 1)) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
% 杨炎龙查看doppler波形 %
% [X,Y] = meshgrid(velocityCoor,distanceCoor);
% figure;surf(Y',X',dopplerSum');title('单个天线的DopplerFFT波形');
% figure(4);imagesc(velocityCoor,distanceCoor,dopplerSum);title('单个天线的DopplerFFT波形');
% set(gca,'YDir','normal');
% xlabel('速度(m/s)');ylabel('距离(m)');zlabel('频谱幅值');

% 结束%
%%% CFAR %%%
guardLen = param.guardLen_range;
winLen = param.winLen_range;
thresholdScale = param.thresholdScale_range;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut = zeros(numADCSamples, dopplerBin_num);

for n = 1: numADCSamples

    CFARDopplerDomainOut(n, :) = cfar_ca(dopplerSum(n, :), dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);

end
guardLen = param.guardLen_doppler;
winLen = param.winLen_doppler;
thresholdScale = param.thresholdScale_doppler;
cfartype = param.cfartype;
CFARRangeDomainOut = zeros(numADCSamples, dopplerBin_num);

for m = 1: dopplerBin_num

    if ismember(1, CFARDopplerDomainOut(:, m))

        CFARRangeDomainOut(:, m) = cfar_sogo( dopplerSum(:, m), numADCSamples, cfartype, thresholdScale, noiseDivShift, guardLen, winLen);

    end

end

CFAROutTemp = CFARDopplerDomainOut & CFARRangeDomainOut;
[row, col] = find(CFAROutTemp(:, :));
% 画出分别找到的所有峰峰值 -- 杨炎龙
% figure;subplot(3,1,1);imshow(CFAROutTemp',[]);
% subplot(3,1,2);imshow(dopplerSum',[]);
% hold on;
% plot(row,col,'ro');
% 结束

peakValue = dopplerSum(row, col);
% if ~isempty(peakValue)
[ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num);
velocityOut = ((peakGrpingCol - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
distanceOut = (peakGrpingRow * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
CFAROut = [peakGrpingRow.', peakGrpingCol.'];

% subplot(3,1,3);imshow(CFAROutTemp',[]);
% hold on; plot(peakGrpingRow,peakGrpingCol,'ro');
% title('二维峰峰值');xlabel("采样数");ylabel("chirp");
%%% AOA FFT %%%
angleBin_num = param.angleBin_num;
MAX_VEL_ENH_PROCESSING = 0;
FinalResult = zeros(length(CFAROut(:, 1)), 5);

if TX_num > 1

    [azimuthOut, elevOut, pcStrc] = AOA2_v1_4(dopplerOut, CFAROut, distanceCoor,velocityCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frame_index, lidarDataFrame,lineId);
%     azimuthOutSize = size(azimuthOut)
%     elevOutSize = size(elevOut)
%     figure;
%     for i=1:azimuthOutSize(1)
%         plot(abs(azimuthOut(i,:)));
%         hold on
%     end
%     subplot(1,4,4);
%     plot(distanceCoor,mmwavedata(1,:,1,1));
    for m = 1: length(CFAROut(:, 1))

        singleAzimuthOut = azimuthOut(m, :);

        [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, angleBin_num);
%         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
        range = (CFAROut(m, 1) * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
        velo = ((CFAROut(m, 2) - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
        xyzInfo = xyzEstimation(maxIdx, range, azimuthOut(m, :), elevOut(m, :), TX_num, angleBin_num);
        FinalResult(m, :) = [range, velo, xyzInfo];

    end

end
    
% else
%     
%     distanceOut = 0;
%     velocityOut = 0;
%     FinalResult = 0;
% 
% end

end