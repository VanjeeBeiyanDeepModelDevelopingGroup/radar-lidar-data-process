%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults_urrsrrNormal function                                 %%%
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

% function [ distanceCoor_vel1, velocityCoor_vel1,velocityCoor_vel2, distanceOut_vel, velocityOut_vel, CFAROut_vel, mmwavedata_vel1, dopplerSum_vel1, dopplerSum_vel2 ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data, frameIndex)
% function [ distanceCoor_vel1, velocityCoor_vel1,velocityCoor_vel2,distanceCoor,velocityCoor, distanceOut_vel, velocityOut_vel, CFAROut_vel, mmwavedata_vel1, dopplerSum, dopplerSum_vel1, dopplerSum_vel2,pcStrc1 ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data,lidarDataFrame, frameIndex,lineId,Temp1_accumulate)
% function [  ] = mmwaveResults_urrsrrNormal_applied(radarDataAll, vel1data, vel2data,lidarDataFrame, frameIndex,lineId,Temp1_accumulate)
function [  ] = mmwaveResults_urrsrrNormal_applied(radarDataAll, vel1data, vel2data,frameIndex)
numADCSamples = param.numADCSamples;
numADCSamples_vel = param.numADCSamples_vel;
numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
RX_num = param.RX_num;
TX_num = param.TX_num;

frame_index = str2num(frameIndex);

[data_num,~,~,~] = size(radarDataAll);
dopplerBin_num = data_num/RX_num;
% ret = saveDataCube(radarDataAll,frame_index);
%%% range FFT %%%
% original part %
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
% 所有数据构成128x512x3xn，即32x4个chirp，512个采样点，3个发送天线，n帧
% 沿着第二个维度，做一维fft
rangeOut = fftshift(fft(rangeOut, numADCSamples, 2),2);

c = param.c;
digOutSampleRate = param.digOutSampleRate;
freqSlopeConst = param.freqSlopeConst;  % 测试那个指令
slope = freqSlopeConst * 1e12;

rangeAbs = abs(rangeOut);

distanceCoor = ((1: numADCSamples) * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
% 返回的是第n帧的第1个天线的第一个chirp的距离FFT
% 更改，返回当前帧的第一个天线的所有chirp的range-FFT
mmwavedata = rangeAbs(:, :, 1, frame_index);

% vel part %
rangeWin_vel = hanning(numADCSamples_vel);
rangeWin_vel = rangeWin_vel(1: (numADCSamples_vel / 2));
RangeWinLen_vel               = length(rangeWin_vel);
RangeWindowCoeffVec_vel       = ones(numADCSamples_vel, 1);
RangeWindowCoeffVec_vel(1:RangeWinLen_vel) = rangeWin_vel;
RangeWindowCoeffVec_vel(numADCSamples_vel-RangeWinLen_vel+1:numADCSamples_vel) = RangeWindowCoeffVec_vel(RangeWinLen_vel:-1:1);
rangeWin_vel   = RangeWindowCoeffVec_vel;

rangeOut_vel1 = bsxfun(@times, vel1data, rangeWin_vel.');
rangeOut_vel1 = fftshift(fft(rangeOut_vel1, numADCSamples_vel, 2),2);
rangeOut_vel2 = bsxfun(@times, vel2data, rangeWin_vel.');
rangeOut_vel2 = fftshift(fft(rangeOut_vel2, numADCSamples_vel, 2),2);

digOutSampleRate_vel = param.digOutSampleRate_vel;
freqSlopeConst_vel = param.freqSlopeConst_vel;  % 测试那个指令
slope_vel = freqSlopeConst_vel * 1e12;

distanceCoor_vel1 = ((1: numADCSamples_vel) * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);

rangeAbs_vel1 = abs(rangeOut_vel1);
mmwavedata_vel1 = rangeAbs_vel1(:, :, frame_index);

%%% doppler FFT %%%
% original part %
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
dopplerOut = fftshift(fft(dopplerIn, dopplerBin_num, 2));
dopplerLog2Abs = 20*log10(abs(dopplerOut));
dopplerSum = sum(dopplerLog2Abs, [3 4]);% 这里应该做beamforming吧，可以做beamforming，可以试试
dopplerSum = squeeze(dopplerSum); 
% 多普勒累积增强
% dopplerSum = SNREnhance(squeeze(sum(abs(dopplerOut),[3,4])),param.USRR_enhance);
% 结束
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

%% 画出doppler图
global nextCnt
% this_frameNum = frame_index+(nextCnt-1)*param.packLength;
% lidarData_frame = lidarDataFrame(:,:,this_frameNum); 
% startNum = 40;
% lidarAngleGrid = (lidarData_frame(1:end-4,1)*256+lidarData_frame(1:end-4,2))/100-90;
% lidarData = lidarData_frame(1:end-4,startNum:end);
% [m,n] = size(lidarData);
% lidarRangeGrid = 0.15*([1:n])-0.4546;
% h=figure(2);figureName = ['帧号：',num2str(this_frameNum)];
% set(h,'name',figureName,'Numbertitle','off')
% subplot(4,2,1);
% figure(3);subplot(3,1,1);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData');hold on
% set(gca,'YDIR','normal');xlim([-30,30]);
% figure(2);subplot(4,2,3);
% imagesc(distanceCoor,velocityCoor,dopplerSum');
% set(gca,'YDIR','normal');hold on
%% 结束
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
dopplerOut_vel1 = fftshift(fft(dopplerIn_vel1, numCirpsPerFrame_vel1, 2));
dopplerLog2Abs_vel1 = abs(dopplerOut_vel1);
dopplerSum_vel1 = sum(dopplerLog2Abs_vel1, 3);
dopplerSum_vel1 = squeeze(dopplerSum_vel1);
% 信噪比增强
% dopplerSum_vel1 = SNREnhance(squeeze(sum(abs(dopplerOut_vel1),[3,4])),param.MRR_enhance);
dopplerIn_vel2 = bsxfun(@times, dopplerIn_vel2, dopplerWin_vel.');
dopplerOut_vel2 = fftshift(fft(dopplerIn_vel2, numCirpsPerFrame_vel2, 2));
dopplerLog2Abs_vel2 = abs(dopplerOut_vel2);
dopplerSum_vel2 = sum(dopplerLog2Abs_vel2, 3);
dopplerSum_vel2 = squeeze(dopplerSum_vel2);
% 信噪比增强
dopplerSum_vel2 = SNREnhance(squeeze(sum(abs(dopplerOut_vel2),[3,4])),param.MRR_enhance);
startFreqConst_vel = param.startFreqConst_vel;
startFreq_vel = startFreqConst_vel * 1e9;
bandwidth_vel = (slope_vel * numADCSamples_vel) / (digOutSampleRate_vel * 1e3);
adcStartTimeConst_vel = param.adcStartTimeConst_vel;
adcStart_vel = adcStartTimeConst_vel * 1e-6;
centerFreq_vel = startFreq_vel + bandwidth_vel * 0.5 + adcStart_vel * slope_vel;
idleTimeConst_vel1 = param.MRR_idleTimeConst_1;
idleTimeConst_vel2 = param.MRR_idleTimeConst_2;
rampEndTime_vel = param.MRR_rampEndTime;
chirpInterval_vel1 = (idleTimeConst_vel1 + rampEndTime_vel) * 1e-6;
chirpInterval_vel2 = (idleTimeConst_vel2 + rampEndTime_vel) * 1e-6;

velocityCoor_vel1 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
velocityCoor_vel2 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel2 * 1000);
% figure(2);
% subplot(4,2,2);
% imagesc(distanceCoor_vel1,velocityCoor_vel1,dopplerSum_vel1');
% set(gca,'YDIR','normal');hold on
% subplot(4,2,4);
% imagesc(distanceCoor_vel1,velocityCoor_vel2,dopplerSum_vel2');
% set(gca,'YDIR','normal');
% 结束%
%%% CFAR %%%
% original part %
guardLen = param.guardLen_doppler;
winLen = param.winLen_doppler;
thresholdScale = param.thresholdScale_doppler;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut = zeros(numADCSamples, dopplerBin_num);

for n = 1: numADCSamples

%     CFARDopplerDomainOut(n, :) = cfar_ca(dopplerSum(n, :), dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
    CFARDopplerDomainOut(n, :) = cfar_ca1D_square( dopplerSum(n,:), winLen, guardLen, thresholdScale, 0);

end
% figure;imagesc(CFARDopplerDomainOut);set(gca,'YDIR', 'normal');
guardLen = param.guardLen_range;
winLen = param.winLen_range;
thresholdScale = param.thresholdScale_range;
cfartype = param.cfartype;
CFARRangeDomainOut = zeros(numADCSamples, dopplerBin_num);

for m = 1: dopplerBin_num

    if ismember(1, CFARDopplerDomainOut(:, m))

%         CFARRangeDomainOut(:, m) = cfar_sogo( dopplerSum(:, m), numADCSamples, cfartype, thresholdScale, noiseDivShift, guardLen, winLen);
        CFARRangeDomainOut(:, m) = cfar_ca1D_square( dopplerSum(:, m), winLen, guardLen, thresholdScale, 0);
    end

end
CFAROutTemp = CFARDopplerDomainOut & CFARRangeDomainOut;
[row, col] = find(CFAROutTemp(:, :));
% figure;imagesc(CFAROutTemp');set(gca,'YDIR', 'normal');
peakValue = dopplerSum([row, col]);
% [ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num);
[ peakGrpingRow, peakGrpingCol ] = cfarFilter0Vel(row, col, dopplerBin_num);
velocityOut = ((peakGrpingCol - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
distanceOut = (peakGrpingRow * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
CFAROut = [peakGrpingRow.', peakGrpingCol.'];
if ~isempty(CFAROut)
    CFAROut = sortrows(CFAROut,1);
end
% CFAROut = [row, col];

% vel CFAR test
MAX_NUM_DET_PER_RANGE_GATE = param.MAX_NUM_DET_PER_RANGE_GATE;
guardLen = param.MRR_guardLen_doppler;
winLen = param.MRR_winLen_doppler;
thresholdScale = param.MRR_thresholdScale_doppler;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);

for n = 1: numADCSamples_vel

%     CFARDopplerDomainOut_vel1(n, :) = cfar_ca(dopplerSum_vel1(n, :), numCirpsPerFrame_vel1, thresholdScale, noiseDivShift, guardLen, winLen);
    CFARDopplerDomainOut_vel1(n, :) = cfar_ca1D_square(dopplerSum_vel1(n, :), winLen, guardLen, thresholdScale, 0);
    
end
% figure;imagesc(CFARDopplerDomainOut_vel1');set(gca,'YDIR', 'normal');
guardLen = param.MRR_guardLen_range;
winLen = param.MRR_winLen_range;
thresholdScale = param.MRR_thresholdScale_range;
cfartype = param.cfartype;
CFARRangeDomainOut_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
% CFARRangeDomainOut_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel2);
distanceOut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
velocityOut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
CFAROut_vel = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);

for m = 1: numCirpsPerFrame_vel1

    if ismember(1, CFARDopplerDomainOut_vel1(:, m))
        
        CFARRangeDomainOut_vel1(:, m) = cfar_ca1D_square( dopplerSum_vel1(:, m), winLen,guardLen , thresholdScale, 0);

    end
    
end


CFAROutTemp_vel1 = CFARDopplerDomainOut_vel1 & CFARRangeDomainOut_vel1;
% figure;imagesc(CFAROutTemp_vel1');set(gca,'YDIR', 'normal');
[row_vel1, col_vel1] = find(CFAROutTemp_vel1(:, :));
peakValue_vel1 = dopplerSum_vel1([row_vel1, col_vel1]);
% [ peakGrpingRow_vel1, peakGrpingCol_vel1 ] = peakPruning(row_vel1, col_vel1, peakValue_vel1, numADCSamples_vel, numCirpsPerFrame_vel1);
[ peakGrpingRow_vel1, peakGrpingCol_vel1 ] = cfarFilter0Vel(row_vel1, col_vel1, numCirpsPerFrame_vel1);
velocityOut_vel = ((peakGrpingCol_vel1 - numCirpsPerFrame_vel1 / 2 - 1) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
distanceOut_vel = (peakGrpingRow_vel1 * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);
CFAROut_vel = [peakGrpingRow_vel1.', peakGrpingCol_vel1.'];
% CFAROut_vel = [row_vel1, col_vel1];
if ~isempty(CFAROut_vel)
    CFAROut_vel = sortrows(CFAROut_vel,1);
end
%%% AOA FFT %%%
angleBin_num = param.angleBin_num;
MAX_VEL_ENH_PROCESSING = 0;
% FinalResult = zeros(length(CFAROut(:, 1)), 5);
%% 计算USRR的角度
fprintf('计算USRR的角度\n');
tic
pcStrc1 = [];
w = linspace(-1,1,angleBin_num); % angle_grid
agl_grid = asin(w)*180/pi; % [-1,1]->[-pi/2,pi/2]
w1 = linspace(-1,1,param.musicBin);
agl_grid_music = asin(w1)*180/pi;
[azimuthOut, elevOut] = AOA2_v1_3(dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING);
% hh=figure;
[numTarget,~] = size(CFAROut);
for m = 1: numTarget

    singleAzimuthOut = azimuthOut(m, :);
%     figure(3);subplot(3,1,2);
%     plot(agl_grid_music,abs(singleAzimuthOut));hold on;
    [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, param.musicBin);
%         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
    range = distanceCoor(CFAROut(m,1));
    vel = velocityCoor(CFAROut(m,2));
    angle = agl_grid_music(maxIdx);
    pcStrc1 = [pcStrc1; [angle,range,vel]];
end
% hold off
% close(hh);
toc
%% 计算MRR的角度
fprintf('计算MRR的角度\n');
tic
[azimuthOut, elevOut] = AOA2_v1_3(dopplerOut_vel1, CFAROut_vel, 1, RX_num, numCirpsPerFrame_vel1, angleBin_num, MAX_VEL_ENH_PROCESSING);

% [rangeAngleMap,dopplerAngleMap] = datacubeProcess(MRR_azimuthDataCube_fft);
[numTarget,~] = size(CFAROut_vel);
% agl_grid = agl_grid-7;
% distanceCoor_vel1 = distanceCoor_vel1+8;
pcStrc2 = [];
for m = 1: numTarget

    singleAzimuthOut = azimuthOut(m, :);
%     figure(3);subplot(3,1,3);
%     plot(agl_grid_music,abs(singleAzimuthOut));hold on;
    [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, param.musicBin);
%         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
    range_vel1 = distanceCoor_vel1(CFAROut_vel(m,1));
    vel_vel1 = velocityCoor_vel1(CFAROut_vel(m,2));
    angle = agl_grid_music(maxIdx);
    pcStrc2 = [pcStrc2; [angle,range_vel1,vel_vel1]];

end
% hold off
toc
% close(hh);
% %% 画出毫米波输出点云
% % USRR图
% if ~isempty(CFAROut)
%     range = distanceCoor(CFAROut(:,1));
%     vel = velocityCoor(CFAROut(:,2));
% %     figure(2);subplot(4,2,3);plot(range,vel,'rp'); % 画在usrr doppler图上
%     
%     % 3d位置
%     angleArr = pcStrc1(:,1);
%     rangeArr = pcStrc1(:,2);
%     pcStrc_usrr = calElevAngLidar(lineId,pcStrc1);
%     y_usrr = -pcStrc_usrr.vertex.x(~isnan(pcStrc_usrr.vertex.x));
%     x_usrr = pcStrc_usrr.vertex.y(~isnan(pcStrc_usrr.vertex.y));
%     z_usrr = pcStrc_usrr.vertex.z(~isnan(pcStrc_usrr.vertex.z));
% %     figure(2);subplot(4,2,1);
% %     figure(3);subplot(3,1,1);plot(angleArr,rangeArr,'rp');
%     
% %     xlim([-30,30]);
% %     xlim([0,140]);
% end
% % MRR图  
% if ~isempty(CFAROut_vel)
%     range_vel1 = distanceCoor_vel1(CFAROut_vel(:,1));
%     vel_vel1 = velocityCoor_vel1(CFAROut_vel(:,2));
% %     figure(2);subplot(4,2,2);plot(range_vel1,vel_vel1,'gp'); % 画在vel1 doppler图上
%     
%     % 3d位置
%     angleArr_mrr = pcStrc2(:,1);
%     rangeArr_mrr = pcStrc2(:,2);
%     pcStrc_mrr = calElevAngLidar(lineId,pcStrc2);
%     y_mrr = -pcStrc_mrr.vertex.x(~isnan(pcStrc_mrr.vertex.x));
%     x_mrr = pcStrc_mrr.vertex.y(~isnan(pcStrc_mrr.vertex.y));
%     z_mrr = pcStrc_mrr.vertex.z(~isnan(pcStrc_mrr.vertex.z));
% %     figure(2);subplot(4,2,1);
% %     figure(3);subplot(3,1,1);plot(angleArr_mrr,rangeArr_mrr,'gp');
% %     xlim([-30,30]);
%     
% %     xlim([0,140]);
% end
%% 计算usrr data cube 3d-fft结果
fprintf('计算usrr data cube 3d-fft结果\n');% 历时 3.258826 秒
tic
% 不做cfar的fft角度谱分析
% [ dataCube3dFFT, ~ ] = AOA2_v1_6( dopplerOut, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING);
% 做cfar的music角度谱分析，不构造矩阵，结果非相干叠加
% [RAMap] = AOA2_v1_7(dopplerOut, CFAROut, TX_num, RX_num, numADCSamples, dopplerBin_num, (TX_num-1)*RX_num, Temp1_accumulate, MAX_VEL_ENH_PROCESSING);
% 不做cfar直接角度谱分析，构造矩阵计算music
% [ RAMap_fromMat_all ] = AOA2_v1_8(dopplerOut, TX_num, RX_num, numADCSamples, dopplerBin_num, (TX_num-1)*RX_num, MAX_VEL_ENH_PROCESSING);
% [carte_RAMap_fromMat,xMat_grid_usrr,yMat_grid_usrr] = polar2carte(RAMap_fromMat_all,agl_grid_music,distanceCoor,1);
% figure;imagesc(xMat_grid_usrr,yMat_grid_usrr,carte_RAMap_fromMat);set(gca,'YDIR','normal');title('usrr range-angle map');
% 做cfar再角分析，构造矩阵计算music
[ RAMap_fromMat_cfar ] = AOA2_v1_9( dopplerOut, CFAROut, TX_num, RX_num, numADCSamples, dopplerBin_num);
% [carte_RAMap_fromMat,xMat_grid_usrr,yMat_grid_usrr] = polar2carte(RAMap_fromMat_cfar,agl_grid_music,distanceCoor,1);
% figure;subplot(1,2,2);imagesc(xMat_grid_usrr,yMat_grid_usrr,carte_RAMap_fromMat);set(gca,'YDIR','normal');title('usrr range-angle map');
% subplot(1,2,1);imagesc(velocityCoor,distanceCoor,dopplerSum);title('usrr range-doppler map');set(gca,'YDIR','normal');
% toc
% [ RAMap_fft, elevOut ] = AOA2_v1_10( dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING);
% 数据累积
% [rangeDoppler_sum,rangeAngle_sum,DopplerAngle_sum] = dataBinning(dataCube3dFFT);
% 画图
% figure(2);subplot(4,2,5);
% imagesc(agl_grid_music,distanceCoor,rangeAngle_sum);hold on
% imagesc(agl_grid_music,distanceCoor,RAMap');hold on
% 把图变换到直角坐标系下
% [carte_RAMap_srr_fft,xgrid_srr_fft,ygrid_srr_fft] = polar2carte(rangeAngle_sum,agl_grid,distanceCoor,1);
% [carte_RAMap_srr_cfarfft,xgrid_srr_cfarfft,ygrid_srr_cfarfft] = polar2carte(RAMap_fft,agl_grid,distanceCoor,1);
[carte_RAMap,xgrid_usrr,ygrid_usrr] = polar2carte(RAMap_fromMat_cfar,agl_grid_music,distanceCoor,1);
% [carte_RAMap,xgrid_usrr,ygrid_usrr] = polar2carte(RAMap',agl_grid,distanceCoor,1);
% imagesc(xgrid_usrr,ygrid_usrr,carte_RAMap);hold on
% if ~isempty(CFAROut)
% %     subplot(4,2,5);plot(angleArr,rangeArr,'rp');
%     subplot(4,2,5);plot(rangeArr.*sin(angleArr*pi/180),rangeArr.*cos(angleArr*pi/180),'rp');
% end
% set(gca,'YDIR','normal');title('usrr range-angle map');
% subplot(4,2,7);
% imagesc(agl_grid,velocityCoor,DopplerAngle_sum');
% set(gca,'YDIR','normal');title('usrr vel-angle map');
% figure(3);subplot(2,1,1);imagesc(distanceCoor,velocityCoor,rangeDoppler_sum');set(gca,'YDIR','normal');
toc
%% 计算mrr data cube 3d-fft结果
fprintf('计算mrr data cube 3d-fft结果\n');% 
tic
% [ dataCube3dFFT_mrr, ~ ] = AOA2_v1_6( dopplerOut_vel1, 1, RX_num, 256, numCirpsPerFrame_vel1, angleBin_num, MAX_VEL_ENH_PROCESSING);
% [RAMap_mrr] = AOA2_v1_7(dopplerOut_vel1, CFAROut_vel, 1, RX_num, 256, numCirpsPerFrame_vel1, angleBin_num, Temp1_accumulate, MAX_VEL_ENH_PROCESSING);
[ RAMap_fromMat_cfar_mrr ] = AOA2_v1_9( dopplerOut_vel1, CFAROut_vel, 1, RX_num, 256, numCirpsPerFrame_vel1);
% [ RAMap_fft_mrr, ~ ] = AOA2_v1_10( dopplerOut_vel1, CFAROut_vel, 1, RX_num, 256, angleBin_num, MAX_VEL_ENH_PROCESSING);
% 数据累积
% [mrr_rangeDoppler_sum,mrr_rangeAngle_sum,mrr_DopplerAngle_sum] = dataBinning(dataCube3dFFT_mrr);
% 画图
% figure(2);subplot(4,2,6);
% imagesc(agl_grid_music,distanceCoor_vel1,RAMap_mrr');hold on
% 把图变换到直角坐标系下
% [carte_RAMap_mrr_fft,xgrid_mrr_fft,ygrid_mrr_fft] = polar2carte(mrr_rangeAngle_sum,agl_grid,distanceCoor_vel1,1);
% [carte_RAMap_mrr_cfarfft,xgrid_mrr_cfarfft,ygrid_mrr_cfarfft] = polar2carte(RAMap_fft_mrr,agl_grid,distanceCoor_vel1,1);
[carte_RAMap_mrr,xgrid_mrr,ygrid_mrr] = polar2carte(RAMap_fromMat_cfar_mrr,agl_grid_music,distanceCoor_vel1,1);
% [carte_RAMap_mrr,xgrid_mrr,ygrid_mrr] = polar2carte(RAMap_mrr',agl_grid,distanceCoor_vel1,1);
% imagesc(xgrid_mrr,ygrid_mrr,carte_RAMap_mrr);hold on
% if ~isempty(CFAROut_vel)
% %     subplot(4,2,6);plot(angleArr_mrr,rangeArr_mrr,'gp');
%     subplot(4,2,6);plot(rangeArr_mrr.*sin(angleArr_mrr*pi/180),rangeArr_mrr.*cos(angleArr_mrr*pi/180),'gp');
% end
% set(gca,'YDIR','normal');title('mrr range-angle map');
% subplot(4,2,8);
% imagesc(agl_grid,velocityCoor_vel1,mrr_DopplerAngle_sum');
% set(gca,'YDIR','normal');title('mrr vel-angle map');
% figure(3);subplot(2,1,2);imagesc(distanceCoor_vel1,velocityCoor_vel1,mrr_rangeDoppler_sum');set(gca,'YDIR','normal');
toc
%% 计算激光雷达光流并统一输入到距离——速度维度
% [lidarFlow] = geneLidarRangeVelMap(lidarData', rangeDoppler_sum, lidarRangeGrid,lidarAngleGrid, distanceCoor, velocityCoor);
% [lidarFlow] = geneLidarRangeVelMap(lidarData', mrr_rangeDoppler_sum, lidarRangeGrid,lidarAngleGrid, distanceCoor_vel1, velocityCoor_vel1);
%% 计算激光雷达的测距点云（使用脉宽修正方法精确测距）lidarData_frame
% fprintf('画激光点云到毫米波雷达RA图上\n');
% tic
% [pcStrc,pcPolar] = lidarRangeMeas(lidarData_frame',lineId+2);
% x = pcStrc.vertex.x(~isnan(pcStrc.vertex.x));
% y = pcStrc.vertex.y(~isnan(pcStrc.vertex.y));
% z = pcStrc.vertex.z(~isnan(pcStrc.vertex.z));
% figure(h);subplot(4,2,1);
% scatter3(x,y,z,4,'b','filled');view([0,0,1]);hold on
% if ~isempty(CFAROut)
%     scatter3(x_usrr,y_usrr,z_usrr,4,'r','filled');view([0,0,1]);
% end
% if ~isempty(CFAROut_vel)
%     scatter3(x_mrr,y_mrr,z_mrr,4,'g','filled');view([0,0,1]);
% end
% subplot(4,2,1);hold off
% xlim([-15,15]);ylim([0,100]);zlim([-2,2]);
% 把激光点云画到毫米波RA图上去
% subplot(4,2,5);scatter(pcPolar(:,1)-90,pcPolar(:,2),3,'r','filled');
% subplot(4,2,6);scatter(pcPolar(:,1)-90,pcPolar(:,2),3,'r','filled');
% subplot(4,2,5);scatter(x,y,3,'y','filled');
% subplot(4,2,6);scatter(x,y,3,'y','filled');
% pause(0.001);
toc
%% music分析结果
fprintf('画出music分析结果\n');
tic
% [rangeSpectrum,angleSpectrum] = spectrumAnalysis_0917(dataCube3dFFT);
h=figure(2);%figureName = ['帧号：',num2str(this_frameNum)];
% subplot(2,2,1);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData');
subplot(2,3,1);plot(distanceCoor,20*log10(sum(dopplerSum',1)));
subplot(2,3,4);plot(distanceCoor_vel1,20*log10(sum(dopplerSum_vel1',1)));

subplot(2,3,2);imagesc(distanceCoor_vel1,velocityCoor,dopplerSum');
subplot(2,3,5);imagesc(distanceCoor_vel1,velocityCoor,dopplerSum_vel1');
set(gca,'YDIR','normal');
% subplot(3,2,1);
subplot(2,3,3);
% [X_usrr,Y_usrr] = meshgrid(xgrid_usrr,ygrid_usrr);
imagesc(xgrid_usrr,ygrid_usrr,carte_RAMap);hold on
% scatter(x,y,3,'r','filled');
% mesh(X_usrr,Y_usrr,carte_RAMap);
set(gca,'YDIR','normal');title('CFAR MUSIC srr range-angle map');
% xlim([5*min(x),5*max(x)]);ylim([0,max(y)]);
% subplot(3,2,2);
subplot(2,3,6);
% [X_mrr,Y_mrr] = meshgrid(xgrid_mrr,ygrid_mrr);
imagesc(xgrid_mrr,ygrid_mrr,carte_RAMap_mrr);hold on
% scatter(x,y,3,'r','filled');
set(gca,'YDIR','normal');title('CFAR MUSIC mrr range-angle map');
% xlim([5*min(x),5*max(x)]);ylim([0,max(y)]);
toc
vars = {'adcStart', 'adcStart_vel', 'adcStartTimeConst', 'adcStartTimeConst_vel', 'agl_grid', 'agl_grid_music', 'angle', 'angleBin_num', 'azimuthOut', 'bandwidth', 'bandwidth_vel', 'c', 'carte_RAMap', 'carte_RAMap_mrr', 'centerFreq', 'centerFreq_vel', 'CFARDopplerDomainOut', 'CFARDopplerDomainOut_vel1', 'CFAROut', 'CFAROut_vel', 'CFAROutTemp', 'CFAROutTemp_vel1', 'CFARRangeDomainOut', 'CFARRangeDomainOut_vel1', 'cfartype', 'chirpInterval', 'chirpInterval_vel1', 'chirpInterval_vel2', 'col', 'col_vel1', 'data_num', 'digOutSampleRate', 'digOutSampleRate_vel', 'distanceCoor', 'distanceCoor_vel1','distanceOut', 'distanceOut_vel', 'dopplerBin_num', 'dopplerIn', 'dopplerIn_vel1', 'dopplerIn_vel2', 'dopplerLog2Abs', 'dopplerLog2Abs_vel1', 'dopplerLog2Abs_vel2', 'dopplerOut', 'dopplerOut_vel1', 'dopplerOut_vel2', 'dopplerSum', 'dopplerSum_vel1', 'dopplerSum_vel2', 'dopplerWin', 'dopplerWin_vel', 'dopplerWindowCoeffVec', 'dopplerWindowCoeffVec_vel', 'dopplerWinLen', 'dopplerWinLen_vel', 'elevOut', 'frame_index', 'frameIndex', 'freqSlopeConst', 'freqSlopeConst_vel', 'guardLen', 'i', 'idleTimeConst', 'idleTimeConst_vel1', 'idleTimeConst_vel2', 'lidarAngleGrid', 'lidarData', 'lidarData_frame', 'lidarDataFrame', 'lidarRangeGrid', 'lineId', 'm', 'MAX_NUM_DET_PER_RANGE_GATE', 'MAX_VEL_ENH_PROCESSING', 'maxIdx', 'mmwavedata', 'mmwavedata_vel1', 'n', 'nextCnt', 'noiseDivShift', 'numADCSamples', 'numADCSamples_vel', 'numCirpsPerFrame_vel1', 'numCirpsPerFrame_vel2', 'numTarget', 'pcPolar', 'pcStrc', 'pcStrc1', 'pcStrc2', 'peakGrpingCol', 'peakGrpingCol_vel1', 'peakGrpingRow', 'peakGrpingRow_vel1', 'peakValue', 'peakValue_vel1', 'radarDataAll', 'RAMap_fromMat_cfar', 'RAMap_fromMat_cfar_mrr', 'rampEndTime', 'rampEndTime_vel', 'range', 'range_vel1', 'rangeAbs', 'rangeAbs_vel1', 'rangeOut', 'rangeOut_vel1', 'rangeOut_vel2', 'rangeWin', 'rangeWin_vel', 'RangeWindowCoeffVec', 'RangeWindowCoeffVec_vel', 'RangeWinLen', 'RangeWinLen_vel', 'row', 'row_vel1', 'RX_num', 'singleAzimuthOut', 'slope', 'slope_vel', 'startFreq', 'startFreq_vel', 'startFreqConst', 'startFreqConst_vel', 'startNum', 'this_frameNum', 'thresholdScale', 'TX_num', 'vel', 'vel1data', 'vel2data', 'vel_vel1', 'velocityCoor', 'velocityCoor_vel1', 'velocityCoor_vel2', 'velocityOut', 'velocityOut_vel', 'w', 'w1', 'winLen', 'x', 'xgrid_mrr', 'xgrid_usrr', 'y', 'ygrid_mrr', 'ygrid_usrr', 'z'};
clear(vars{:})
end