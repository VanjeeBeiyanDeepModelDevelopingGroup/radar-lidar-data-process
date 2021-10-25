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

function [ azimuthOut, elevOut ] = AOA2_v1_3( dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING)
    if isempty(CFAROut)
        azimuthOut = [];
        elevOut = [];
        return;
    end
    azimuthOut = zeros(length(CFAROut(:, 1)), angleBin_num);
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
%              save Temp1_.mat Temp1
             azimuthOut(m,:) = fftshift(fft(Temp1));
             elevOut(m,:) = fftshift(fft(Temp2));
% %% music �㷨
% %              [P1,f1] = pmusic(Temp1,30,64,4); % input,p,nfft,fs,nwin,overlap
%              [P1,f1] = pmusic(Temp1,[Inf,1.1],64,64,3); % Window length = 7
%              azimuthOut(m, :) = P1;
% %              [P2,f2] = pmusic(Temp2,30,64,4); % input,p,nfft,fs,nwin,overlap
%              [P2,f2] = pmusic(Temp1,[Inf,1.1],64,64,3); % Window length = 7
%              elevOut(m, :) = P2;
%              figure(5);
%              subplot(2,2,1);plot(20*log10(abs(azimuthFFTResult)));
%              subplot(2,2,3);plot(f1,P1);legend('fft','music');title('usrr_�ź��ӿռ�ά��Ϊ30/64');
%              a=1;
         else

             Temp1 = zeros(1, angleBin_num);
             Temp1(1: TX_num * RX_num) = reshape(dopplerOut(CFAROut(m, 1), CFAROut(m, 2), :), [1, TX_num * RX_num]);
             Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
             Temp1 = dopplerCompensation(Temp1, CFAROut(m, 2), TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);

             azimuthOut(m,:) = fftshift(fft(Temp1));
% %% music �㷨
% %              [P1,f1] = pmusic(Temp1,15,64,4); % input,p,nfft,fs,nwin,overlap
%              [P1,f1] = pmusic(Temp1,[Inf,1.1],64,64,3); % Window length = 7
% 
%              azimuthOut(m, :) = P1;
%              subplot(2,2,2);plot(20*log10(abs(azimuthFFTResult)));
%              subplot(2,2,4);plot(f1,P1);legend('fft','music');title('mrr_�ź��ӿռ�ά��Ϊ4/64');
%              a=1;
         end

     end
     

end