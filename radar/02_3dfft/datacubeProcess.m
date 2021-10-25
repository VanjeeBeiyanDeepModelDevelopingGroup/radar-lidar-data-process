function [RAmap,DAmap] = dataceubeProces(radarDataAll, vel1data, vel2data, lidarDataFrame_singL, num2str(i),lineId)
%% 基本参数
USRR_dopplerBin = 32;
USRR_RXNum = 3;
[USRR_rxMultDopplerBin, USRR_rangeBin, USRR_TXNum, USRR_frameNum] = size(radarDataAll);
[MRR_rxMultDopplerBin, MRR_rangeBin, MRR_frameNum] = size(vel1data);
%% 将数据打包成rangeBin*dopplerBin*antennaBin形式的datacube
USRR_dataCube = zeros(numADCSamples, dopplerBin_num, RX_num, TX_num);% 512x32x4x3
for n = 1: USRR_rangeBin
    for i = 1: USRR_TXNum
        for m = 1: USRR_RXNum
            USRR_dataCube(n, :, m, i) = rangeOut((0: (dopplerBin_num - 1)) * RX_num + m, n, i, frame_index);
        end
    end
end


end