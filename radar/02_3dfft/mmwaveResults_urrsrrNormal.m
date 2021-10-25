%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults_urrsrrNormal function                                 %%%
%%%     radarDataAll- ���ײ��״�����                                        %%%
%%%     vel1data, vel2data- fast chirp�״����ݣ�slow chirp�״�����          %%%
%%%     lidarDataFrame- �����״�����                                        %%%
%%%     frameIndex- ��鿴��֡����λ��                                      %%%
%%%                                                                         %%%
%%%     distanceCoor- ���������̶�                                        %%%
%%%     velocityCoor- �������ٶȿ̶�                                        %%%
%%%     distanceOut- ѡ��֡��CFAR��ɸ�ĵ�ľ�����Ϣ                         %%%
%%%     velocityOut- ѡ��֡��CFAR��ɸ�ĵ���ٶ���Ϣ                         %%%
%%%     FinalResult- ѡ��֡��CFAR��ɸ�ĵ�ľ��룬�ٶȣ�xyz��������Ϣ        %%%
%%%     mmwavedata- һάFFT����ľ���ֵ                                     %%%
%%%     dopplerSum- doppler fft����ļӺ�ֵ                                 %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.05.12 version 1.3                            %%%
%%%     �޸ģ����Ӷ��ַ������ٶȵĴ���                                      %%%
%%%           ����ѡ����methodSign��־λ��������                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [ distanceCoor_vel1, velocityCoor_vel1,velocityCoor_vel2, distanceOut_vel, velocityOut_vel, CFAROut_vel, mmwavedata_vel1, dopplerSum_vel1, dopplerSum_vel2 ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data, frameIndex)
function [ distanceCoor_vel1, velocityCoor_vel1,velocityCoor_vel2,distanceCoor,velocityCoor, distanceOut_vel, velocityOut_vel, CFAROut_vel, mmwavedata_vel1, dopplerSum, dopplerSum_vel1, dopplerSum_vel2,pcStrc1 ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data,lidarDataFrame, frameIndex,lineId)
%% ���ļ��е�ȫ�ֱ���
global CFAROut
global CFAROut_vel


numADCSamples = param.numADCSamples;
numADCSamples_vel = param.numADCSamples_vel;
numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
RX_num = param.RX_num;
TX_num = param.TX_num;

% frame_index = str2num(frameIndex);
frame_index = 1;
[data_num,~,~,~] = size(radarDataAll);
dopplerBin_num = data_num/RX_num;

% ret = saveDataCube(radarDataAll,frame_index);
%%% ȥ��ֱ����������ֵ���֣� %%%
%%% ��Ч�����ÿ�ɾ�� ������ 20210927 %%%
% mean_chirp = mean(radarDataAll,2);
% radarDataAll_ = radarDataAll - mean_chirp;
% radarDataAll = radarDataAll_;
% mean_chirp_vel1 = mean(vel1data,2);
% vel1data_ = vel1data - mean_chirp_vel1;
% vel1data = vel1data_;
% mean_chirp_vel2 = mean(vel2data,2);
% vel2data_ = vel2data - mean_chirp_vel2;
% vel2data = vel2data_;
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
% ��������512ά
rangeOut = bsxfun(@times, radarDataAll, rangeWin.');
% rangeOut_size = size(rangeOut)
% �������ݹ���128x512x3xn����32x4��chirp��512�������㣬3���������ߣ�n֡
% ���ŵڶ���ά�ȣ���һάfft
rangeOut = fftshift(fft(rangeOut, numADCSamples, 2),2);
% rangeOut = fft(rangeOut, numADCSamples, 2);% �ȸɵ�shift�������Ͼ�����һά�Ȳ���Ҫshift���ɣ���

c = param.c;
digOutSampleRate = param.digOutSampleRate;
freqSlopeConst = param.freqSlopeConst;  % �����Ǹ�ָ��
slope = freqSlopeConst * 1e12;

rangeAbs = abs(rangeOut);

distanceCoor = ((1: numADCSamples) * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
% ���ص��ǵ�n֡�ĵ�1�����ߵĵ�һ��chirp�ľ���FFT
% ���ģ����ص�ǰ֡�ĵ�һ�����ߵ�����chirp��range-FFT
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
% rangeOut_vel1 = fft(rangeOut_vel1, numADCSamples_vel, 2);% �ȸɵ�shift�������Ͼ�����һά�Ȳ���Ҫshift���ɣ���
rangeOut_vel2 = bsxfun(@times, vel2data, rangeWin_vel.');
rangeOut_vel2 = fftshift(fft(rangeOut_vel2, numADCSamples_vel, 2),2);
% rangeOut_vel2 = fft(rangeOut_vel2, numADCSamples_vel, 2);% �ȸɵ�shift�������Ͼ�����һά�Ȳ���Ҫshift���ɣ���

digOutSampleRate_vel = param.digOutSampleRate_vel;
freqSlopeConst_vel = param.freqSlopeConst_vel;  % �����Ǹ�ָ��
slope_vel = freqSlopeConst_vel * 1e12;

distanceCoor_vel1 = ((1: numADCSamples_vel) * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);

rangeAbs_vel1 = abs(rangeOut_vel1);
mmwavedata_vel1 = rangeAbs_vel1(:, :, frame_index);

%%% doppler FFT %%%
% original part %
dopplerWin = hanning(dopplerBin_num);  % 32��chirp
dopplerWin = dopplerWin(1: (dopplerBin_num / 2));
dopplerWinLen               = length(dopplerWin);
dopplerWindowCoeffVec       = ones(dopplerBin_num, 1);
dopplerWindowCoeffVec(1:dopplerWinLen) = dopplerWin;
dopplerWindowCoeffVec(dopplerBin_num-dopplerWinLen+1:dopplerBin_num) = dopplerWindowCoeffVec(dopplerWinLen:-1:1);
dopplerWin   = dopplerWindowCoeffVec;
% urr ��128*512*3������תΪ512*32*4*3�ı�׼��ʽ %
dopplerIn = zeros(numADCSamples, dopplerBin_num, RX_num, TX_num);% 512x32x4x3
for n = 1: numADCSamples
    for i = 1: TX_num
        for m = 1: RX_num
            dopplerIn(n, :, m, i) = rangeOut((0: (dopplerBin_num - 1)) * RX_num + m, n, i, frame_index);
        end
    end
end
dopplerIn = bsxfun(@times, dopplerIn, dopplerWin.');% �Ӵ�
% dopplerOut = fftshift(fft(dopplerIn, dopplerBin_num, 2));
%%% ��ԭ����FFT֮��ֱ��shift��Ϊ���ɵ�0Ƶ�ʸ�����ֵ %%%
%%% Ч�����ÿ�ɾ���������� 20210928 %%%
dopplerOut = fft(dopplerIn, dopplerBin_num, 2);
dopplerOut(:,[1,2,32],:,:) = 0;
dopplerOut = fftshift(dopplerOut);

dopplerLog2Abs = 20*log10(abs(dopplerOut));
dopplerSum = sum(dopplerLog2Abs, [3 4]);% ����Ӧ����beamforming�ɣ�������beamforming����������
dopplerSum = squeeze(dopplerSum); 
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

%% ����dopplerͼ
global nextCnt
this_frameNum = frame_index+(nextCnt-1)*param.packLength;
lidarData_frame = lidarDataFrame(:,:,this_frameNum); 
startNum = 40;
lidarAngleGrid = (lidarData_frame(1:end-4,1)*256+lidarData_frame(1:end-4,2))/100-90;
lidarData = lidarData_frame(1:end-4,startNum:end);
[m,n] = size(lidarData);
lidarRangeGrid = 0.15*([1:n]);
h=figure(2);figureName = ['֡�ţ�',frameIndex];
set(h,'name',figureName,'Numbertitle','off')
% subplot(4,2,1);
figure(3);subplot(3,1,1);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData');hold on
set(gca,'YDIR','normal');xlim([-30,30]);
figure(2);subplot(4,2,3);
h2_6=imagesc(distanceCoor,velocityCoor,dopplerSum');
% set(h2_6,'WindowButtonDownFcn',{@ButttonDownFcn});
title('usrr vel-range map');
set(gca,'YDIR','normal');hold on
%% ����
% vel part %
dopplerWin_vel = hanning(numCirpsPerFrame_vel1);
dopplerWin_vel = dopplerWin_vel(1: (numCirpsPerFrame_vel1 / 2));
dopplerWinLen_vel               = length(dopplerWin_vel);
dopplerWindowCoeffVec_vel       = ones(numCirpsPerFrame_vel1, 1);
dopplerWindowCoeffVec_vel(1:dopplerWinLen_vel) = dopplerWin_vel;
dopplerWindowCoeffVec_vel(numCirpsPerFrame_vel1-dopplerWinLen_vel+1:numCirpsPerFrame_vel1) = dopplerWindowCoeffVec_vel(dopplerWinLen_vel:-1:1);
dopplerWin_vel   = dopplerWindowCoeffVec_vel;
% mrr ��128*512*3������תΪ256*128*4�ı�׼��ʽ %
dopplerIn_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1, RX_num);
dopplerIn_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1, RX_num);
for m = 1: RX_num
    dopplerIn_vel1(:, :, m) = rangeOut_vel1((0: (numCirpsPerFrame_vel1 - 1)) * RX_num + m, :, frame_index).';
    dopplerIn_vel2(:, :, m) = rangeOut_vel2((0: (numCirpsPerFrame_vel1 - 1)) * RX_num + m, :, frame_index).';
end

dopplerIn_vel1 = bsxfun(@times, dopplerIn_vel1, dopplerWin_vel.');
% dopplerOut_vel1 = fftshift(fft(dopplerIn_vel1, numCirpsPerFrame_vel1, 2));
%%% ��ԭ����FFT֮��ֱ��shift��Ϊ���ɵ�0Ƶ�ʸ�����ֵ %%%
%%% Ч�����ÿ�ɾ���������� 20210928 %%%
dopplerOut_vel1 = fft(dopplerIn_vel1, numCirpsPerFrame_vel1, 2);
dopplerOut_vel1(:,[1,2,128],:) = 0;
dopplerOut_vel1 = fftshift(dopplerOut_vel1);

dopplerLog2Abs_vel1 = abs(dopplerOut_vel1);
dopplerSum_vel1 = 20*log10(sum(dopplerLog2Abs_vel1, 3));
dopplerSum_vel1 = squeeze(dopplerSum_vel1);

dopplerIn_vel2 = bsxfun(@times, dopplerIn_vel2, dopplerWin_vel.');
% dopplerOut_vel2 = fftshift(fft(dopplerIn_vel2, numCirpsPerFrame_vel2, 2));
%%% ��ԭ����FFT֮��ֱ��shift��Ϊ���ɵ�0Ƶ�ʸ�����ֵ %%%
%%% Ч�����ÿ�ɾ���������� 20210928 %%%
dopplerOut_vel2 = fft(dopplerIn_vel2, numCirpsPerFrame_vel2, 2);
dopplerOut_vel2(:,[1,2,128],:) = 0;
dopplerOut_vel2 = fftshift(dopplerOut_vel2);

dopplerLog2Abs_vel2 = abs(dopplerOut_vel2);
dopplerSum_vel2 = 20*log10(sum(dopplerLog2Abs_vel2, 3));
dopplerSum_vel2 = squeeze(dopplerSum_vel2);

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
figure(2);
subplot(4,2,2);
imagesc(distanceCoor_vel1,velocityCoor_vel1,dopplerSum_vel1');
set(gca,'YDIR','normal');hold on
title('mrr vel-range map');
subplot(4,2,4);
imagesc(distanceCoor_vel1,velocityCoor_vel2,dopplerSum_vel2');
set(gca,'YDIR','normal');
title('mrr vel-range map');
% ����%
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
[ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num);
velocityOut = ((peakGrpingCol - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
distanceOut = (peakGrpingRow * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
CFAROut = [peakGrpingRow.', peakGrpingCol.'];
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
[ peakGrpingRow_vel1, peakGrpingCol_vel1 ] = peakPruning(row_vel1, col_vel1, peakValue_vel1, numADCSamples_vel, numCirpsPerFrame_vel1);
velocityOut_vel = ((peakGrpingCol_vel1 - numCirpsPerFrame_vel1 / 2 - 1) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
distanceOut_vel = (peakGrpingRow_vel1 * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);
CFAROut_vel = [peakGrpingRow_vel1.', peakGrpingCol_vel1.'];
% CFAROut_vel = [row_vel1, col_vel1];
%%% AOA FFT %%%
angleBin_num = param.angleBin_num;
MAX_VEL_ENH_PROCESSING = 0;
% FinalResult = zeros(length(CFAROut(:, 1)), 5);
%% ����USRR�ĽǶ�

pcStrc1 = [];
w = linspace(-1,1,angleBin_num); % angle_grid
agl_grid = asin(w)*180/pi+7; % [-1,1]->[-pi/2,pi/2]
[azimuthOut, elevOut] = AOA2_v1_3(dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING);
% hh=figure;
[numTarget,~] = size(CFAROut);
for m = 1: numTarget

    singleAzimuthOut = azimuthOut(m, :);
    figure(3);subplot(3,1,2);
    plot(agl_grid,abs(singleAzimuthOut),'o-');
    hold on
    [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, angleBin_num);
%         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
    range = distanceCoor(CFAROut(m,1));
    vel = velocityCoor(CFAROut(m,2));
    angle = agl_grid(maxIdx);
    pcStrc1 = [pcStrc1; [angle,range,vel]];
end
figure(3);subplot(3,1,2);hold off
% close(hh);
    
%% ����MRR�ĽǶ�
[azimuthOut, elevOut] = AOA2_v1_3(dopplerOut_vel1, CFAROut_vel, 1, RX_num, numCirpsPerFrame_vel1, angleBin_num, MAX_VEL_ENH_PROCESSING);
% [rangeAngleMap,dopplerAngleMap] = datacubeProcess(MRR_azimuthDataCube_fft);
[numTarget,~] = size(CFAROut_vel);
% agl_grid = agl_grid-7;
% distanceCoor_vel1 = distanceCoor_vel1+8;
pcStrc2 = [];
for m = 1: numTarget

    singleAzimuthOut = azimuthOut(m, :);
    figure(3);subplot(3,1,3);
    plot(agl_grid,abs(singleAzimuthOut),'o-');
    hold on
    [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, angleBin_num);
%         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
    range_vel1 = distanceCoor_vel1(CFAROut_vel(m,1));
    vel_vel1 = velocityCoor_vel1(CFAROut_vel(m,2));
    angle = agl_grid(maxIdx);
    pcStrc2 = [pcStrc2; [angle,range_vel1,vel_vel1]];

end
figure(3);subplot(3,1,3);hold off
% close(hh);
%% �������ײ��������
% USRRͼ
if ~isempty(CFAROut)
    range = distanceCoor(CFAROut(:,1));
    vel = velocityCoor(CFAROut(:,2));
    figure(2);subplot(4,2,3);plot(range,vel,'rp'); % ����usrr dopplerͼ��
    title('usrr vel-range map');
    hold off
    % 3dλ��
    angleArr = pcStrc1(:,1);
    rangeArr = pcStrc1(:,2);
    pcStrc_usrr = calElevAngLidar(lineId,pcStrc1);
    y_usrr = -pcStrc_usrr.vertex.x(~isnan(pcStrc_usrr.vertex.x));
    x_usrr = pcStrc_usrr.vertex.y(~isnan(pcStrc_usrr.vertex.y));
    z_usrr = pcStrc_usrr.vertex.z(~isnan(pcStrc_usrr.vertex.z));
    figure(3);subplot(3,1,1);plot(angleArr,rangeArr,'rp');
%     for kk=1:length(rangeArr)
%         figure(3);subplot(3,1,1);plot(agl_grid,rangeArr(kk)*ones(1,length(agl_grid)),'r-');
%     end
%     xlim([-30,30]);
%     xlim([0,140]);
end
% MRRͼ  
if ~isempty(CFAROut_vel)
    range_vel1 = distanceCoor_vel1(CFAROut_vel(:,1));
    vel_vel1 = velocityCoor_vel1(CFAROut_vel(:,2));
    figure(2);subplot(4,2,2);plot(range_vel1,vel_vel1,'gp'); % ����vel1 dopplerͼ��
    title('mrr vel-range map');
    hold off
    % 3dλ��
    angleArr_mrr = pcStrc2(:,1);
    rangeArr_mrr = pcStrc2(:,2);
    pcStrc_mrr = calElevAngLidar(lineId,pcStrc2);
    y_mrr = -pcStrc_mrr.vertex.x(~isnan(pcStrc_mrr.vertex.x));
    x_mrr = pcStrc_mrr.vertex.y(~isnan(pcStrc_mrr.vertex.y));
    z_mrr = pcStrc_mrr.vertex.z(~isnan(pcStrc_mrr.vertex.z));
%     figure(2);subplot(4,2,1);
    figure(3);subplot(3,1,1);plot(angleArr_mrr,rangeArr_mrr,'gp');
%     for kk=1:length(rangeArr_mrr)
%         figure(3);subplot(3,1,1);plot(agl_grid,rangeArr_mrr(kk)*ones(1,length(agl_grid)),'g');
%     end
    xlim([-30,30]);
    hold off
%     xlim([0,140]);
end
%% ����usrr data cube 3d-fft���

[ dataCube3dFFT, ~ ] = AOA2_v1_6( dopplerOut, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING);
% �����ۻ�
[rangeDoppler_sum,rangeAngle_sum,DopplerAngle_sum] = dataBinning(dataCube3dFFT);
% ��ͼ
figure(2);subplot(4,2,5);
imagesc(agl_grid,distanceCoor,rangeAngle_sum);hold on
if ~isempty(CFAROut)
    subplot(4,2,5);plot(angleArr,rangeArr,'rp');hold on
end
set(gca,'YDIR','normal');title('usrr range-angle map');
xlim([-30,30]);
subplot(4,2,7);
imagesc(agl_grid,velocityCoor,DopplerAngle_sum');
set(gca,'YDIR','normal');title('usrr vel-angle map');

% figure(3);subplot(2,1,1);imagesc(distanceCoor,velocityCoor,rangeDoppler_sum');set(gca,'YDIR','normal');
%% ����mrr data cube 3d-fft���
tic
[ dataCube3dFFT_mrr, ~ ] = AOA2_v1_6( dopplerOut_vel1, 1, RX_num, 256, numCirpsPerFrame_vel1, angleBin_num, MAX_VEL_ENH_PROCESSING);
toc
% �����ۻ�
[mrr_rangeDoppler_sum,mrr_rangeAngle_sum,mrr_DopplerAngle_sum] = dataBinning(dataCube3dFFT_mrr);
% ��ͼ
figure(2);subplot(4,2,6);
imagesc(agl_grid,distanceCoor_vel1,mrr_rangeAngle_sum);hold on

if ~isempty(CFAROut_vel)
    subplot(4,2,6);plot(angleArr_mrr,rangeArr_mrr,'gp');hold on
end
set(gca,'YDIR','normal');title('mrr range-angle map');
xlim([-30,30]);ylim([0,100]);
subplot(4,2,8);
imagesc(agl_grid,velocityCoor_vel1,mrr_DopplerAngle_sum');
set(gca,'YDIR','normal');title('mrr vel-angle map');
% figure(3);subplot(2,1,2);imagesc(distanceCoor_vel1,velocityCoor_vel1,mrr_rangeDoppler_sum');set(gca,'YDIR','normal');

%% ���㼤���״������ͳһ���뵽���롪���ٶ�ά��
% [lidarFlow] = geneLidarRangeVelMap(lidarData', rangeDoppler_sum, lidarRangeGrid,lidarAngleGrid, distanceCoor, velocityCoor);
% [lidarFlow] = geneLidarRangeVelMap(lidarData', mrr_rangeDoppler_sum, lidarRangeGrid,lidarAngleGrid, distanceCoor_vel1, velocityCoor_vel1);
%% ���㼤���״�Ĳ����ƣ�ʹ����������������ȷ��ࣩlidarData_frame
[pcStrc,pcPolar] = lidarRangeMeas(lidarData_frame',lineId+2);
x = pcStrc.vertex.x(~isnan(pcStrc.vertex.x));
y = pcStrc.vertex.y(~isnan(pcStrc.vertex.y));
z = pcStrc.vertex.z(~isnan(pcStrc.vertex.z));
figure(h);subplot(4,2,1);
scatter3(x,y,z,4,'b','filled');view([0,0,1]);hold on
if ~isempty(CFAROut)
    scatter3(x_usrr,y_usrr,z_usrr,4,'r','filled');view([0,0,1]);
end
if ~isempty(CFAROut_vel)
    scatter3(x_mrr,y_mrr,z_mrr,4,'g','filled');view([0,0,1]);
end
title('���ƶԱ�');
xlim([-15,15]);ylim([0,100]);zlim([-2,2]);hold off
% �����ƻ������ײ��Ƕ�ͼ��ȥ
subplot(4,2,5);scatter(pcPolar(:,1)-90,pcPolar(:,2),3,'y','filled');hold off
subplot(4,2,6);scatter(pcPolar(:,1)-90,pcPolar(:,2),3,'r','filled');hold off
% title(point_cloud,'����','FontWeight','bold','FontSize',15)

% %% ��ӻص����������ѡ��
% select_points = datacursormode(h);
% set(select_points, 'UpdateFcn', {@myUpdateFcn});

%% ��ͣ1ms
pause(0.001);
figure(3);hold off
figure(h);hold off
end

% function ButttonDownFcn(hObject, ~)
% global CFAROut
% global CFAROut_vel
% pt = get(hObject,'CurrentPoint');
% x = pt(1,1);
% y = pt(1,2);
% % fprintf('x=%f,y=%f\n',x,y);
% % subplot(4,2,1)
% hold on;plot(x,y,'ko','MarkerSize',5);
% end
% 
% function pts = myUpdateFcn(source, event_obj)
% %% �������ѡ��Ļص�����
% % TODO: �ں��ٸ�
% pos=source.Position;
% 
% childrens=event_obj.Target.Parent.Parent.Children;
% 
% if event_obj.Target.Parent==childrens(7)
%     a=1;
%     hold(childrens(5),'on')
%     plot(childrens(5),1:100,2:2:200)
% end
% 
% end