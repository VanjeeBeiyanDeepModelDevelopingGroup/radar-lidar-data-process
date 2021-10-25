%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults_CRT function                                          %%%
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

function [ distanceCoor_vel1, velocityCoor_vel1, distanceOut_vel, velocityOut_vel, mmwavedata_vel1, dopplerSum_vel1 ] = mmwaveResults_CRT(radarDataAll, vel1data, vel2data, frameIndex)
% function [ distanceCoor, velocityCoor, distanceOut, velocityOut, FinalResult, mmwavedata, dopplerSum,  ] = mmwaveResults_v1_2(radarDataAll, lidarDataFrame, frameIndex)

numADCSamples = param.numADCSamples;
% numChirpsPerFrame = param.numChirpsPerFrame;
numADCSamples_vel = param.numADCSamples_vel;
numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
RX_num = param.RX_num;
TX_num = param.TX_num;
% numDataPerFrame = numChirpsPerFrame * RX_num * numADCSamples * 2;

frame_index = str2num(frameIndex);

% data_num = 128
[data_num,~,~,~] = size(radarDataAll);
dopplerBin_num = data_num/RX_num;
% data_num

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
% �������ݹ���128x512x3x8����32x4��chirp��512�������㣬3���������ߣ�8֡
% ���ŵڶ���ά�ȣ���һάfft
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
rangeOut_vel1 = fft(rangeOut_vel1, numADCSamples_vel, 2);
rangeOut_vel2 = bsxfun(@times, vel2data, rangeWin_vel.');
rangeOut_vel2 = fft(rangeOut_vel2, numADCSamples_vel, 2);

digOutSampleRate_vel = param.digOutSampleRate_vel;
freqSlopeConst_vel = param.freqSlopeConst_vel;  % �����Ǹ�ָ��
slope_vel = freqSlopeConst_vel * 1e12;

distanceCoor_vel1 = ((1: numADCSamples_vel) * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);

rangeAbs_vel1 = abs(rangeOut_vel1);
mmwavedata_vel1 = rangeAbs_vel1(:, :, frame_index);

% mmwavedataSize = size(mmwavedata)
% ����һ�²��
% guardLen = 2;
% winLen = 16;
% thresholdScale = 40;
% noiseDivShift = ceil(log2(2 * winLen));
% cfarResult = cfar_ca( mmwavedata(1,:), numADCSamples, thresholdScale, noiseDivShift, guardLen, winLen);
% rangeFFTwaveCluster = splitRangeFFT(mmwavedata(1,:),cfarResult);

%%% doppler FFT %%%
% original part %
dopplerWin = hanning(dopplerBin_num);  % 32��chirp
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
% ��������angle-fft֮ǰ��xcube
% figure(2);
% dopplerLog2Abs_3d(:,:,1:4) = dopplerLog2Abs(:,:,:,1);
% dopplerLog2Abs_3d(:,:,5:8) = dopplerLog2Abs(:,:,:,2);
% dopplerLog2Abs_3d(:,:,9:12) = dopplerLog2Abs(:,:,:,3);
% for i=1:12
%     subplot(3,4,i);
%     imagesc(dopplerLog2Abs_3d(:,:,i));set(gca,'YDir','normal');title(num2str(i));
% end
% dopplerLog2AbsSize = size(dopplerLog2Abs)
dopplerSum = sum(dopplerLog2Abs, [3 4]);% ����Ӧ����beamforming�ɣ�������beamforming����������
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
% �������鿴doppler���� %
% [X,Y] = meshgrid(velocityCoor,distanceCoor);
% figure;surf(Y',X',dopplerSum');title('�������ߵ�DopplerFFT����');
% figure(4);imagesc(velocityCoor,distanceCoor,dopplerSum);title('�������ߵ�DopplerFFT����');
% set(gca,'YDir','normal');
% xlabel('�ٶ�(m/s)');ylabel('����(m)');zlabel('Ƶ�׷�ֵ');

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
idleTimeConst_vel2 = param.MRR_idleTimeConst_2;
rampEndTime_vel = param.MRR_rampEndTime;
chirpInterval_vel1 = (idleTimeConst_vel1 + rampEndTime_vel) * 1e-6;
chirpInterval_vel2 = (idleTimeConst_vel2 + rampEndTime_vel) * 1e-6;

velocityCoor_vel1 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
velocityCoor_vel2 = (((- numCirpsPerFrame_vel1 / 2): (numCirpsPerFrame_vel1 / 2 - 1)) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel2 * 1000);


% ����%
%%% CFAR %%%
% original part %
guardLen = param.guardLen_doppler;
winLen = param.winLen_doppler;
thresholdScale = param.thresholdScale_doppler;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut = zeros(numADCSamples, dopplerBin_num);

for n = 1: numADCSamples

    CFARDopplerDomainOut(n, :) = cfar_ca(dopplerSum(n, :), dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);

end
guardLen = param.guardLen_range;
winLen = param.winLen_range;
thresholdScale = param.thresholdScale_range;
cfartype = param.cfartype;
CFARRangeDomainOut = zeros(numADCSamples, dopplerBin_num);

for m = 1: dopplerBin_num

    if ismember(1, CFARDopplerDomainOut(:, m))

        CFARRangeDomainOut(:, m) = cfar_sogo( dopplerSum(:, m), numADCSamples, cfartype, thresholdScale, noiseDivShift, guardLen, winLen);

    end

end

CFAROutTemp = CFARDopplerDomainOut & CFARRangeDomainOut;
[row, col] = find(CFAROutTemp(:, :));


% �����ֱ��ҵ������з��ֵ -- ������
% figure;subplot(3,1,1);imshow(CFAROutTemp',[]);
% subplot(3,1,2);imshow(dopplerSum',[]);
% hold on;
% plot(row,col,'ro');
% ����

peakValue = dopplerSum(row, col);
[ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num);
velocityOut = ((peakGrpingCol - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
distanceOut = (peakGrpingRow * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
CFAROut = [peakGrpingRow.', peakGrpingCol.'];

% subplot(3,1,3);imshow(CFAROutTemp',[]);
% hold on; plot(peakGrpingRow,peakGrpingCol,'ro');
% title('��ά���ֵ');xlabel("������");ylabel("chirp");

% vel part %
MAX_NUM_DET_PER_RANGE_GATE = param.MAX_NUM_DET_PER_RANGE_GATE;
guardLen = param.guardLen_doppler;
winLen = param.winLen_doppler;
thresholdScale = param.thresholdScale_doppler;
noiseDivShift = ceil(log2(2 * winLen));
CFARDopplerDomainOut_vel1 = zeros(numADCSamples_vel, numCirpsPerFrame_vel1);
% CFARDopplerDomainOut_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel2);
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
% CFARRangeDomainOut_vel2 = zeros(numADCSamples_vel, numCirpsPerFrame_vel2);
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
    
%     if ismember(1, CFARDopplerDomainOut_vel2(:, m))
% 
%         CFARRangeDomainOut_vel2(:, m) = cfar_sogo( dopplerSum_vel2(:, m), numADCSamples_vel, cfartype, thresholdScale, noiseDivShift, guardLen, winLen);
% 
%     end

end


CFAROutTemp_vel1 = CFARDopplerDomainOut_vel1 & CFARRangeDomainOut_vel1;
% dopplerSum_vel1 = squeeze(sum(abs(dopplerIn_vel1),3));

% [row_vel1, col_vel1] = find(CFAROut_vel(:, :));
% peakValue_vel1 = dopplerSum_vel1(row_vel1, col_vel1);
% [ peakGrpingRow_vel1, peakGrpingCol_vel1 ] = peakPruning(row_vel1, col_vel1, peakValue_vel1, numADCSamples_vel, numCirpsPerFrame_vel1);
% velocityOut_vel = ((peakGrpingCol_vel1 - numCirpsPerFrame_vel1 / 2 - 1) * c * 3600) / (2 * numCirpsPerFrame_vel1 * centerFreq_vel * chirpInterval_vel1 * 1000);
% distanceOut_vel = (peakGrpingRow_vel1 * c * digOutSampleRate_vel * 1e3) / (2 * slope_vel * numADCSamples_vel);
% CFAROut_vel = [peakGrpingRow_vel1.', peakGrpingCol_vel1.'];

%%% Draw Graphs %%%


% %%% AOA FFT %%%
% angleBin_num = param.angleBin_num;
% MAX_VEL_ENH_PROCESSING = 0;
% FinalResult = zeros(length(CFAROut(:, 1)), 5);
% 
% if TX_num > 1
%     
%     [azimuthOut, elevOut] = AOA2_v1_4(dopplerOut, CFAROut, distanceCoor,velocityCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frame_index, lidarDataFrame);
% %     azimuthOutSize = size(azimuthOut)
% %     elevOutSize = size(elevOut)
% %     figure;
% %     for i=1:azimuthOutSize(1)
% %         plot(abs(azimuthOut(i,:)));
% %         hold on
% %     end
% %     subplot(1,4,4);
% %     plot(distanceCoor,mmwavedata(1,:,1,1));
%     for m = 1: length(CFAROut(:, 1))
% 
%         singleAzimuthOut = azimuthOut(m, :);
% 
%         [ maxIdx, ~ ] = powerAndMax(singleAzimuthOut, angleBin_num);
% %         [maxIdx,~] = cfar_ca(singleAzimuthOut, dopplerBin_num, thresholdScale, noiseDivShift, guardLen, winLen);
%         range = (CFAROut(m, 1) * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
%         velo = ((CFAROut(m, 2) - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
%         xyzInfo = xyzEstimation(maxIdx, range, azimuthOut(m, :), elevOut(m, :), TX_num, angleBin_num);
%         FinalResult(m, :) = [range, velo, xyzInfo];
% 
%     end
%     
% end

end