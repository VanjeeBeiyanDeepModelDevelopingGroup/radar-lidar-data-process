%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults_CRT function                                          %%%
%%%     radarDataAll- 毫米波雷达数据                                        %%%
%%%     vel1data, vel2data- fast chirp雷达数据，slow chirp雷达数据          %%%
%%%     lidarDataFrame- 激光雷达数据                                        %%%
%%%     frameIndex- 想查看的帧序列位置                                      %%%
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

function [ distanceCoor_vel1, velocityCoor_vel1, distanceOut_vel, velocityOut_vel, mmwavedata_vel1, dopplerSum_vel1 ] = mmwaveResults_CRT_mrr(vel1data, vel2data, frameIndex)

numADCSamples_vel = param.numADCSamples_vel;
numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
RX_num = param.RX_num;

frame_index = str2num(frameIndex);

%%% range FFT %%%
% vel part %
rangeWin_vel = hanning(numADCSamples_vel);
rangeWin_vel = rangeWin_vel(1: (numADCSamples_vel / 2));
RangeWinLen_vel               = length(rangeWin_vel);
RangeWindowCoeffVec_vel       = ones(numADCSamples_vel, 1);
RangeWindowCoeffVec_vel(1:RangeWinLen_vel) = rangeWin_vel;
RangeWindowCoeffVec_vel(numADCSamples_vel-RangeWinLen_vel+1:numADCSamples_vel) = RangeWindowCoeffVec_vel(RangeWinLen_vel:-1:1);
rangeWin_vel   = RangeWindowCoeffVec_vel;

rangeOut_vel1 = bsxfun(@times, vel1data, rangeWin_vel.');
rangeOut_vel1 = fft(rangeOut_vel1, numADCSamples_vel, 2);
rangeOut_vel2 = bsxfun(@times, vel2data, rangeWin_vel.');
rangeOut_vel2 = fft(rangeOut_vel2, numADCSamples_vel, 2);

c = param.c;
digOutSampleRate_vel = param.digOutSampleRate_vel;
freqSlopeConst_vel = param.freqSlopeConst_vel;  % 测试那个指令
slope_vel = freqSlopeConst_vel * 1e12;

distanceCoor_vel1 = ((1: numADCSamples_vel) * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);

rangeAbs_vel1 = abs(rangeOut_vel1);
mmwavedata_vel1 = rangeAbs_vel1(:, :, frame_index);

%%% doppler FFT %%%
% vel part %
dopplerWin_vel = hanning(numCirpsPerFrame_vel1);
dopplerWin_vel = dopplerWin_vel(1: (numCirpsPerFrame_vel1 / 2));
dopplerWinLen_vel               = length(dopplerWin_vel);
dopplerWindowCoeffVec_vel       = ones(numCirpsPerFrame_vel1, 1);
dopplerWindowCoeffVec_vel(1:dopplerWinLen_vel) = dopplerWin_vel;
dopplerWindowCoeffVec_vel(numCirpsPerFrame_vel1-dopplerWinLen_vel+1:numCirpsPerFrame_vel1) = dopplerWindowCoeffVec_vel(dopplerWinLen_vel:-1:1);
dopplerWin_vel   = dopplerWindowCoeffVec_vel;

dopplerIn_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1, RX_num);
dopplerIn_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1, RX_num);

for m = 1: RX_num

    dopplerIn_vel1(:, :, m) = rangeOut_vel1((0: (numCirpsPerFrame_vel1 - 1)) * RX_num + m, :, frame_index).';
    dopplerIn_vel2(:, :, m) = rangeOut_vel2((0: (numCirpsPerFrame_vel1 - 1)) * RX_num + m, :, frame_index).';

end

dopplerIn_vel1 = bsxfun(@times, dopplerIn_vel1, dopplerWin_vel.');
dopplerOut_vel1 = fftshift(fft(dopplerIn_vel1, numCirpsPerFrame_vel1, 2), 2);
dopplerLog2Abs_vel1 = log2(abs(dopplerOut_vel1));
dopplerSum_vel1 = sum(dopplerLog2Abs_vel1, 3);
dopplerSum_vel1 = squeeze(dopplerSum_vel1);

dopplerIn_vel2 = bsxfun(@times, dopplerIn_vel2, dopplerWin_vel.');
dopplerOut_vel2 = fftshift(fft(dopplerIn_vel2, numCirpsPerFrame_vel2, 2), 2);
dopplerLog2Abs_vel2 = log2(abs(dopplerOut_vel2));
dopplerSum_vel2 = sum(dopplerLog2Abs_vel2, 3);
dopplerSum_vel2 = squeeze(dopplerSum_vel2);

startFreqConst_vel = param.startFreqConst_vel;
startFreq_vel = startFreqConst_vel * 1e9;
bandwidth_vel = (slope_vel * numADCSamples_vel) / (digOutSampleRate_vel * 1e3);
adcStartTimeConst_vel = param.adcStartTimeConst_vel;
adcStart_vel = adcStartTimeConst_vel * 1e-6;
centerFreq_vel = startFreq_vel + bandwidth_vel * 0.5 + adcStart_vel * slope_vel;
idleTimeConst_vel1 = param.MRR_idleTimeConst_1;
% idleTimeConst_vel2 = param.MRR_idleTimeConst_2;
rampEndTime_vel = param.MRR_rampEndTime;
chirpInterval_vel1 = (idleTimeConst_vel1 + rampEndTime_vel) * 1e-6;
% chirpInterval_vel2 = (idleTimeConst_vel2 + rampEndTime_vel) * 1e-6;

velocityCoor_vel1 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
% velocityCoor_vel2 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel2 * 1000);

% CFAR %
% vel part %
MAX_NUM_DET_PER_RANGE_GATE = param.MAX_NUM_DET_PER_RANGE_GATE;
guardLen = param.guardLen_doppler;
winLen = param.winLen_doppler;
thresholdScale = param.thresholdScale_doppler;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
CFARDopplerDomainOut_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel2);
DetObj_num = 0;
disambVelIndx = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);

for n = 1: numADCSamples_vel

    CFARDopplerDomainOut_vel1(n, :) = cfar_ca(dopplerSum_vel1(n, :), numCirpsPerFrame_vel1, thresholdScale, noiseDivShift, guardLen, winLen);
%     CFARDopplerDomainOut_vel2(n, :) = cfar_ca(dopplerSum_vel2(n, :), numCirpsPerFrame_vel2, thresholdScale, noiseDivShift, guardLen, winLen);
    if ismember(1, CFARDopplerDomainOut_vel1(n, :))
        [ CFARDopplerDomainOut_vel1(n, :), peaknum, disambVelIndx(n, :) ] = MaxVel(CFARDopplerDomainOut_vel1(n, :), dopplerSum_vel1(n, :), dopplerSum_vel2(n, :), numCirpsPerFrame_vel1, MAX_NUM_DET_PER_RANGE_GATE);
        DetObj_num = DetObj_num + peaknum;
    end
    
end

guardLen = param.guardLen_range;
winLen = param.winLen_range;
thresholdScale = param.thresholdScale_range;
cfartype = param.cfartype;
CFARRangeDomainOut_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
distanceOut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
velocityOut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
CFAROut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
numDet = 0;

for m = 1: numCirpsPerFrame_vel1

    if ismember(1, CFARDopplerDomainOut_vel1(:, m))
        
        CFARRangeDomainOut_vel1(:, m) = cfar_sogo( dopplerSum_vel1(:, m), numADCSamples_vel, cfartype, thresholdScale, noiseDivShift, guardLen, winLen);
        if ismember(1, CFARRangeDomainOut_vel1(:, m))
            [ CFARRangeDomainOut_vel1(:, m), numDet, distanceOut_vel(:, m), velocityOut_vel(:, m), CFAROut_vel(:, m) ] = MaxVel_2D(CFARRangeDomainOut_vel1(:, m), CFARDopplerDomainOut_vel1(:, m), dopplerSum_vel1(:, m), disambVelIndx(:, m), numADCSamples_vel, numDet, m);
        end
            
    end
    

end

CFAROutTemp_vel1 = CFARDopplerDomainOut_vel1 & CFARRangeDomainOut_vel1;

%%% Draw Graphs %%%

end