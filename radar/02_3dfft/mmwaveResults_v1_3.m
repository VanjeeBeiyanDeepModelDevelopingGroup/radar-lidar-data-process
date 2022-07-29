%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     mmwaveResults function                                              %%%
%%%     radarDataAll- ���ײ��״�����                                        %%%
%%%     lidarDataFrame- �����״�����                                        %%%
%%%     frameIndex- ��鿴��֡����λ��                                      %%%
%%%     methodSign- ������õ���������ٶȵı�ʶλ                          %%%
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
% mmwavedataSize = size(mmwavedata)
% ����һ�²��
% guardLen = 2;
% winLen = 16;
% thresholdScale = 40;
% noiseDivShift = ceil(log2(2 * winLen));
% cfarResult = cfar_ca( mmwavedata(1,:), numADCSamples, thresholdScale, noiseDivShift, guardLen, winLen);
% rangeFFTwaveCluster = splitRangeFFT(mmwavedata(1,:),cfarResult);

%%% doppler FFT %%%
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

% ����%
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
% �����ֱ��ҵ������з��ֵ -- ������
% figure;subplot(3,1,1);imshow(CFAROutTemp',[]);
% subplot(3,1,2);imshow(dopplerSum',[]);
% hold on;
% plot(row,col,'ro');
% ����

peakValue = dopplerSum(row, col);
% if ~isempty(peakValue)
[ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num);
velocityOut = ((peakGrpingCol - dopplerBin_num / 2 - 1) * c * 3600) / (2 * dopplerBin_num * TX_num * centerFreq * chirpInterval * 1000);
distanceOut = (peakGrpingRow * c * digOutSampleRate * 1e3) / (2 * slope * numADCSamples);
CFAROut = [peakGrpingRow.', peakGrpingCol.'];

% subplot(3,1,3);imshow(CFAROutTemp',[]);
% hold on; plot(peakGrpingRow,peakGrpingCol,'ro');
% title('��ά���ֵ');xlabel("������");ylabel("chirp");
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