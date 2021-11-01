function [azimuthOut_rodnet_music] = AOA2_v1_7(dopplerOut, CFAROut, TX_num, RX_num, rangeBin_num, dopplerBin_num, angleBin_num, Temp1_accumulate, MAX_VEL_ENH_PROCESSING)
%AOA2_v1_7 ʹ��rodnet�ķ�������RA MAP
%   ����radar cube��RDMap cfar�����������cfar�Ľ���ĽǶ�Ƶ��ƴ��һ��ͼ��
L=360;
azimuthOut_rodnet_music = zeros(L,rangeBin_num);
azimuthOut_rodnet_fft = zeros(angleBin_num,rangeBin_num);
% �ۻ���֡�Ķ������ź�
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
     % ע�����Ӧ�ù���������Ϊchannel�ź�
     Temp1_accumulate = [Temp1_accumulate,Temp1'];
     [~,cols] = size(Temp1_accumulate);
%      if cols > musicTimeWindow
%          % pop ����һ��
%          Temp1_accumulate(:,1) = [];
%          spectOut_music_frameAccumulate = musicAlg(Temp1_accumulate,L);
%      end
     spectOut_Music = musicAlg(Temp1',L);
     spectOut_fft = fftshift(fft((Temp1)));
     range = CFAROut(i,1);
     % ע���������������dopplerά�ȵ��׷��������ɵ���
     % ��ʵ��Ӧ����dopplerά����ֱ�ӹ���doppler*angle��matrix������matrix����music�㷨���з������
     azimuthOut_rodnet_music(:,range) = azimuthOut_rodnet_music(:,range)+spectOut_Music';
     azimuthOut_rodnet_fft(:,range) = azimuthOut_rodnet_fft(:,range)+20*log10(abs(spectOut_fft'));
end

end