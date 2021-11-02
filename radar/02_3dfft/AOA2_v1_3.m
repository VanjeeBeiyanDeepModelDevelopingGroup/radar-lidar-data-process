%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AOA algorithm                                                                            %%%
%%%     dopplerOut- doppler FFT�Ľ������Ϊ����������                                            %%%  
%%%                 ����ά�ȷֱ�Ϊ numADCSamples, dopplerBin_num, RX_num, TX_num, frame_num      %%% 
%%%     CFAROut- CFAR��õ���ֵ��ÿ֡������                                                      %%%
%%%              һ��֡����Ԫ�����飬ÿ��Ԫ�������еĵ�һ����range�������ڶ���Ϊdoppler����      %%%
%%%     TX_num- ������������                                                                     %%%
%%%     RX_num- ������������                                                                     %%%
%%%     dopplerBin_num- doppler FFT bin����                                                      %%%
%%%     angleBin_num- �Ƕȹ���FFT bin����                                                        %%%
%%%     MAX_VEL_ENH_PROCESSING- �Ƿ�ʹ��MAX VEL���ܵı�־λ                                      %%%
%%%     azimuthOut- ˮƽ�Ƕ�3D FFT���                                                           %%%
%%%     elevOut- ��ֱ�Ƕ�3D FFT���                                                              %%%
%%%                                                                                              %%%
%%%     Created by ��α� 2021.02.19 version 1.3                                                 %%%
%%%     �޸Ĳ��֣� ȥ����֡ѭ�����㣬����ѡ��֡���е�֡����                                      %%%
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