function [azimuthOut_rodnet_music] = AOA2_v1_7(dopplerOut, CFAROut, TX_num, RX_num, rangeBin_num, dopplerBin_num, angleBin_num, Temp1_accumulate, MAX_VEL_ENH_PROCESSING)
%AOA2_v1_7 使用rodnet的方法生成RA MAP
%   输入radar cube和RDMap cfar结果，仅仅将cfar的结果的角度频谱拼在一张图里
L=360;
azimuthOut_rodnet_music = zeros(L,rangeBin_num);
azimuthOut_rodnet_fft = zeros(angleBin_num,rangeBin_num);
% 累积多帧的多天线信号
musicTimeWindow = param.musicTimeWindow;
% Temp1_accumulate = [];
for i =1:length(CFAROut(:,1))
    if TX_num == 3
        Temp1 = zeros(1, angleBin_num);
        Temp2 = zeros(1, angleBin_num);

        Temp3 = reshape(dopplerOut(CFAROut(i, 1), CFAROut(i, 2), :, :), [1, TX_num * RX_num]);
        Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
        Temp3 = dopplerCompensation(Temp3, CFAROut(i, 2), TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
        Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num));
    else
        Temp1 = zeros(1, angleBin_num);
        Temp1(1: TX_num * RX_num) = reshape(dopplerOut(CFAROut(i, 1), CFAROut(i, 2), :), [1, TX_num * RX_num]);
        Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
        Temp1 = dopplerCompensation(Temp1, CFAROut(i, 2), TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
    end
     % 注意这里，应该构造列向量为channel信号
     Temp1_accumulate = [Temp1_accumulate,Temp1'];
     [~,cols] = size(Temp1_accumulate);
%      if cols > musicTimeWindow
%          % pop 掉第一列
%          Temp1_accumulate(:,1) = [];
%          spectOut_music_frameAccumulate = musicAlg(Temp1_accumulate,L);
%      end
     spectOut_Music = musicAlg(Temp1',L);
     spectOut_fft = fftshift(fft((Temp1)));
     range = CFAROut(i,1);
     % 注意这里，等于是做了doppler维度的谱分析结果相干叠加
     % 事实上应该在doppler维度上直接构造doppler*angle的matrix，将该matrix送入music算法进行方向估计
     azimuthOut_rodnet_music(:,range) = azimuthOut_rodnet_music(:,range)+spectOut_Music';
     azimuthOut_rodnet_fft(:,range) = azimuthOut_rodnet_fft(:,range)+20*log10(abs(spectOut_fft'));
end

end