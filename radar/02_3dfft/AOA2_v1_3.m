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

function [ azimuthOut_music, elevOut ] = AOA2_v1_3( dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING)
    if isempty(CFAROut)
        azimuthOut = [];
        azimuthOut_music = [];
        elevOut = [];
        return;
    end
    L = 360;
    azimuthOut = zeros(length(CFAROut(:, 1)), angleBin_num);
    azimuthOut_music = zeros(length(CFAROut(:, 1)), L);
    elevOut = zeros(length(CFAROut(:, 1)), angleBin_num);
     for m = 1: length(CFAROut(:, 1))
         if TX_num == 3

             Temp1 = zeros(1, angleBin_num);
             Temp2 = zeros(1, angleBin_num);

             Temp3 = reshape(dopplerOut(CFAROut(m, 1), CFAROut(m, 2), :, :), [1, TX_num * RX_num]);
             Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
             Temp3 = dopplerCompensation(Temp3, CFAROut(m, 2), TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
             Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num));
             Temp2(1: RX_num) = Temp3((((TX_num - 1) * RX_num) + 1): (TX_num * RX_num));
             azimuthOut(m,:) = fftshift(fft(Temp1));
             azimuthOut_music(m,:) = musicAlg(Temp3(1: ((TX_num - 1) * RX_num))',L);
             elevOut(m,:) = fftshift(fft(Temp2));
         else

             Temp1 = zeros(1, angleBin_num);
             Temp1(1: TX_num * RX_num) = reshape(dopplerOut(CFAROut(m, 1), CFAROut(m, 2), :), [1, TX_num * RX_num]);
             Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
             Temp1 = dopplerCompensation(Temp1, CFAROut(m, 2), TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);

             azimuthOut(m,:) = fftshift(fft(Temp1));
             azimuthOut_music(m,:) = musicAlg(Temp1(1: (TX_num * RX_num))',L);
         end

     end
     

end