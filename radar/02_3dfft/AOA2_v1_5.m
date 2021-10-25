%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AOA algorithm                                                                            %%%
%%%     dopplerOut- doppler FFT�Ľ������Ϊ����������                                            %%%  
%%%                 ����ά�ȷֱ�Ϊ numADCSamples, dopplerBin_num, RX_num, TX_num, frame_num      %%% 
%%%     CFAROut- CFAR��õ���ֵ��ÿ֡������                                                      %%%
%%%              һ��֡����Ԫ�����飬ÿ��Ԫ�������еĵ�һ����range�������ڶ���Ϊdoppler����      %%%
%%%     distanceCoor- ��������                                                                   %%%
%%%     TX_num- ������������                                                                     %%%
%%%     RX_num- ������������                                                                     %%%
%%%     numADCSamples- range bin����                                                             %%%
%%%     dopplerBin_num- doppler FFT bin����                                                      %%%
%%%     angleBin_num- �Ƕȹ���FFT bin����                                                        %%%
%%%     MAX_VEL_ENH_PROCESSING- �Ƿ�ʹ��MAX VEL���ܵı�־λ                                      %%%
%%%     azimuthOut- ˮƽ�Ƕ�3D FFT���                                                           %%%
%%%     elevOut- ��ֱ�Ƕ�3D FFT���                                                              %%%
%%%                                                                                              %%%
%%%     Created by ��α� 2021.03.18 version 1.4                                                 %%%
%%%     �޸Ĳ��֣� ȥ��cfar���ɸѡ�����������ݽ���angle fft                                     %%%
%%%                ��range-angle heat map                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ azimuthOut, elevOut, pcStrc ] = AOA2_v1_5( dopplerOut, CFAROut, distanceCoor,velocityCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frameIndex,lidarDataFrame,lineId)
    azimuthOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
    elevOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
     for r = 1: numADCSamples
         
         for c = 1: dopplerBin_num

             if TX_num == 3

                 Temp1 = zeros(1, angleBin_num);
                 Temp2 = zeros(1, angleBin_num);

                 Temp3 = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp3 = dopplerCompensation(Temp3, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
                 Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num)); % 8��ֵ
                 Temp2(1: RX_num) = Temp3((((TX_num - 1) * RX_num) + 1): (TX_num * RX_num)); % 4��ֵ
                 %% �Ӵ�
%                  angleWin = hanning(angleBin_num);  % 32��chirp
%                  angleWin = angleWin(1: (angleBin_num / 2));
%                  angleWinLen               = length(angleWin);
%                  angleWindowCoeffVec       = ones(angleBin_num, 1);
%                  angleWindowCoeffVec(1:angleWinLen) = angleWin;
%                  angleWindowCoeffVec(angleBin_num-angleWinLen+1:angleBin_num) = angleWindowCoeffVec(angleWinLen:-1:1);
%                  angleWin = angleWindowCoeffVec;
%                  Temp1 = bsxfun(@times, Temp1, angleWin.');
                 azimuthOut(r, :, c) = fftshift(fft(Temp1));
                 elevOut(r, :, c) = fftshift(fft(Temp2));

             else

                 Temp1 = zeros(1, angleBin_num);
                 Temp1(1: TX_num * RX_num) = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp1 = dopplerCompensation(Temp1, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);

                 azimuthOut(r, :, c) = fftshift(fft(Temp1));

             end
         
         end

     end
     
     w = linspace(-1,1,angleBin_num); % angle_grid
     agl_grid = asin(w)*180/pi+5; % [-1,1]->[-pi/2,pi/2]
     Xpow = abs(azimuthOut(:, :, :)).^2;
     %% �ٶȹ���ȡĿ��ز�����01
     Xpow1 = generateHiSNRRAmap(dopplerOut,Xpow);
     if Xpow1(1,1) ~= 0
         Xpow1 = 20*log(Xpow1);
     end
     %% �źŻ������ǿ����02
     [rr,dd,aa] = size(Xpow);
     Xpow_mult = ones([rr,dd]);
     for ang = 1:aa
         Xpow_mult = Xpow_mult.*(20*log10(Xpow(:,:,ang)));
     end
     if Xpow_mult(1,1) ~= 0
         Xpow_mult = 20*log10(Xpow_mult);
     end
     %% ֱ�ӰѺ��ײ��ͼ���ԭʼͼƬ����һ��
     % �����״�����
     lidarData = lidarDataFrame(1:end-1,:,frameIndex);
     lidarAngleGrid = (lidarData(:,1)*256+lidarData(:,2))/100-90;
     startNum = 40;
     lidarData = lidarData(:,startNum:end);  % ȥ��ÿ��ǰ13����
     [~,n] = size(lidarData);
     lidarRangeGrid = 0.15*([1:n]);
	 % ��һ������
     % ���ײ��״�
     radarRawMap = Xpow1;  % ����ȡ��һ���ĺ��ײ��ز�����
     maxRadar = max(radarRawMap(:));
     minRadar = min(radarRawMap(:));
     radarRawMap_norm = (radarRawMap-minRadar)/(maxRadar-minRadar);
     % ���ײ��״�
     maxLidar = max(lidarData(:));
     minLidar = min(lidarData(:));
     lidarData_norm = (lidarData-minLidar)/(maxLidar-minLidar);
     % ��ֵ
     [outRawFusionMap,radarRawInLidarMap,new_angleGrid,new_rangeGrid] = fusionV0_1(radarRawMap_norm,lidarData_norm,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
     %% ת��Ϊ�ѿ�������ϵ
%      [carteXSNR,x,y] = polar2carte(Xsnr, agl_grid, distanceCoor, 2);
     %% ��ʾ
     
end