%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AOA algorithm                                                                            %%%
%%%     dopplerOut- doppler FFT的结果，作为函数的输入                                            %%%  
%%%                 矩阵维度分别为 numADCSamples, dopplerBin_num, RX_num, TX_num, frame_num      %%% 
%%%     CFAROut- CFAR后得到峰值的每帧行列数                                                      %%%
%%%              一共帧数个元胞数组，每个元胞数组中的第一列是range索引，第二列为doppler索引      %%%
%%%     TX_num- 发射天线数量                                                                     %%%
%%%     RX_num- 接收天线数量                                                                     %%%
%%%     dopplerBin_num- doppler FFT bin容量                                                      %%%
%%%     angleBin_num- 角度估计FFT bin容量                                                        %%%
%%%     MAX_VEL_ENH_PROCESSING- 是否使能MAX VEL功能的标志位                                      %%%
%%%     azimuthOut- 水平角度3D FFT结果                                                           %%%
%%%     elevOut- 垂直角度3D FFT结果                                                              %%%
%%%                                                                                              %%%
%%%     Created by 李嘉宝 2021.02.19 version 1.3                                                 %%%
%%%     修改部分： 去除多帧循环运算，仅对选定帧进行单帧运算                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ azimuthOut, elevOut] = AOA2_v1_2( dopplerOut, CFAROut, distanceCoor,velocityCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frameIndex,lidarDataFrame,lineId)
% function [ azimuthOut, elevOut, agl_grid, ] = AOA2_v1_4( dopplerOut, CFAROut, distanceCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frameIndex,lidarDataFrame)
    azimuthOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
    elevOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
%     azimuthOut = zeros(numADCSamples, 512, dopplerBin_num);
%     elevOut = zeros(numADCSamples, 512, dopplerBin_num);
%      figure(12);
     for r = 1: numADCSamples
         
         for c = 1: dopplerBin_num

             if TX_num == 3

                 Temp1 = zeros(1, angleBin_num);
                 Temp2 = zeros(1, angleBin_num);

                 Temp3 = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp3 = dopplerCompensation(Temp3, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
                 Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num)); % 8个值
                 Temp2(1: RX_num) = Temp3((((TX_num - 1) * RX_num) + 1): (TX_num * RX_num)); % 4个值
                 %% 加窗
%                  angleWin = hanning(angleBin_num);  % 32个chirp
%                  angleWin = angleWin(1: (angleBin_num / 2));
%                  angleWinLen               = length(angleWin);
%                  angleWindowCoeffVec       = ones(angleBin_num, 1);
%                  angleWindowCoeffVec(1:angleWinLen) = angleWin;
%                  angleWindowCoeffVec(angleBin_num-angleWinLen+1:angleBin_num) = angleWindowCoeffVec(angleWinLen:-1:1);
%                  angleWin = angleWindowCoeffVec;
%                  Temp1 = bsxfun(@times, Temp1, angleWin.');
                 azimuthOut(r, :, c) = fftshift(fft(Temp1));
%                  [music,f1] = pmusic(Temp1,32,angleBin_num,'whole');
%                  azimuthOut(r, :, c) = 20*log10(abs(music));
%                  plot(abs(azimuthOut(r, :, c)));hold on
%                  azimuthOut(r, :, c) = fft(Temp1);
                 elevOut(r, :, c) = fftshift(fft(Temp2));
%                  [elevOut(r, :, c),f2] = pmusic(Temp2,32,angleBin_num,'whole');

             else

                 Temp1 = zeros(1, angleBin_num);
                 Temp1(1: TX_num * RX_num) = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp1 = dopplerCompensation(Temp1, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);

                 azimuthOut(r, :, c) = fftshift(fft(Temp1));
%                 [azimuthOut(r, :, c),f1] = pmusic(Temp1,32,angleBin_num,'whole');

             end
         
         end

     end

end