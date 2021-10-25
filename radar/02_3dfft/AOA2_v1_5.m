%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AOA algorithm                                                                            %%%
%%%     dopplerOut- doppler FFT的结果，作为函数的输入                                            %%%  
%%%                 矩阵维度分别为 numADCSamples, dopplerBin_num, RX_num, TX_num, frame_num      %%% 
%%%     CFAROut- CFAR后得到峰值的每帧行列数                                                      %%%
%%%              一共帧数个元胞数组，每个元胞数组中的第一列是range索引，第二列为doppler索引      %%%
%%%     distanceCoor- 距离坐标                                                                   %%%
%%%     TX_num- 发射天线数量                                                                     %%%
%%%     RX_num- 接收天线数量                                                                     %%%
%%%     numADCSamples- range bin数量                                                             %%%
%%%     dopplerBin_num- doppler FFT bin容量                                                      %%%
%%%     angleBin_num- 角度估计FFT bin容量                                                        %%%
%%%     MAX_VEL_ENH_PROCESSING- 是否使能MAX VEL功能的标志位                                      %%%
%%%     azimuthOut- 水平角度3D FFT结果                                                           %%%
%%%     elevOut- 垂直角度3D FFT结果                                                              %%%
%%%                                                                                              %%%
%%%     Created by 李嘉宝 2021.03.18 version 1.4                                                 %%%
%%%     修改部分： 去除cfar结果筛选，对整个数据进行angle fft                                     %%%
%%%                画range-angle heat map                                                        %%%
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
     %% 速度过滤取目标回波——01
     Xpow1 = generateHiSNRRAmap(dopplerOut,Xpow);
     if Xpow1(1,1) ~= 0
         Xpow1 = 20*log(Xpow1);
     end
     %% 信号互相关增强——02
     [rr,dd,aa] = size(Xpow);
     Xpow_mult = ones([rr,dd]);
     for ang = 1:aa
         Xpow_mult = Xpow_mult.*(20*log10(Xpow(:,:,ang)));
     end
     if Xpow_mult(1,1) ~= 0
         Xpow_mult = 20*log10(Xpow_mult);
     end
     %% 直接把毫米波和激光原始图片插在一起
     % 激光雷达数据
     lidarData = lidarDataFrame(1:end-1,:,frameIndex);
     lidarAngleGrid = (lidarData(:,1)*256+lidarData(:,2))/100-90;
     startNum = 40;
     lidarData = lidarData(:,startNum:end);  % 去掉每包前13个数
     [~,n] = size(lidarData);
     lidarRangeGrid = 0.15*([1:n]);
	 % 归一化操作
     % 毫米波雷达
     radarRawMap = Xpow1;  % 可以取不一样的毫米波回波数据
     maxRadar = max(radarRawMap(:));
     minRadar = min(radarRawMap(:));
     radarRawMap_norm = (radarRawMap-minRadar)/(maxRadar-minRadar);
     % 毫米波雷达
     maxLidar = max(lidarData(:));
     minLidar = min(lidarData(:));
     lidarData_norm = (lidarData-minLidar)/(maxLidar-minLidar);
     % 插值
     [outRawFusionMap,radarRawInLidarMap,new_angleGrid,new_rangeGrid] = fusionV0_1(radarRawMap_norm,lidarData_norm,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
     %% 转换为笛卡尔坐标系
%      [carteXSNR,x,y] = polar2carte(Xsnr, agl_grid, distanceCoor, 2);
     %% 显示
     
end