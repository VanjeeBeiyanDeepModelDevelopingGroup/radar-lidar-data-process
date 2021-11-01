% 取range下的slow_time*channel矩阵输入music算法计算angle谱响应
% 特点：遍历所有距离，取doppler-angle层进行music计算 -- 2021-10-30
function [ spectrumAngle ] = AOA2_v1_8( dopplerOut, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING)
    % angleBin_num 输入真实天线个数
    L = 360;
    if TX_num == 3
        angleNum = (TX_num - 1) * RX_num;
    else
        angleNum = TX_num * RX_num;
    end
    dopplerOut3d = zeros(numADCSamples, dopplerBin_num,angleNum);
     for r = 1: numADCSamples         
         for c = 1: dopplerBin_num
             if c >= dopplerBin_num/2 && c <= dopplerBin_num/2+2
                continue;
             end
             if TX_num == 3
%                  dopplerOut3d = zeros(numADCSamples, dopplerBin_num,(TX_num - 1) * RX_num);
                 Temp3 = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp3 = dopplerCompensation(Temp3, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
                 dopplerOut3d(r, c, :) = Temp3(1: ((TX_num - 1) * RX_num)); % 8个值
             else
%                  dopplerOut3d = zeros(numADCSamples, dopplerBin_num, TX_num*RX_num);
                 Temp1 = zeros(1, angleBin_num);
                 Temp1(1: TX_num * RX_num) = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp1 = dopplerCompensation(Temp1, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
                 dopplerOut3d(r,c,:) = Temp1(1:TX_num*RX_num);
             end
         end
     end
     [numADCSamples,~,~] = size(dopplerOut3d);
     spectrumAngle = zeros(numADCSamples,L);
     for r=1:numADCSamples
         signalMat = squeeze(dopplerOut3d(r,:,:));
         spectrumAngle(r,:) = musicAlg(signalMat',L);
     end
end