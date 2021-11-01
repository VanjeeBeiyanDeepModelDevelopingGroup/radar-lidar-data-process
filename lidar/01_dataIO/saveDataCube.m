%
% 存储data cube
%
function [ret] = saveDataCube(radarDataAll,frame_index)
%     [dopplerBin_num,numADCSamples,TX_num, frameNum] = size(radarDataAll);
    RX_num = 4;
    TX_num = 3;
    dopplerBin_num = 32;
    numADCSamples = 512;
    dopplerIn = zeros(numADCSamples,dopplerBin_num,RX_num,TX_num);
    for n = 1: numADCSamples
        for i = 1: TX_num
            for m = 1: RX_num
                dopplerIn(n, :, m, i) = radarDataAll((0: (dopplerBin_num - 1)) * RX_num + m, n, i, frame_index);
            end
        end
    end
    datacube = reshape(dopplerIn,[numADCSamples,dopplerBin_num,RX_num*TX_num]);
    save datacube.mat datacube
    if ~isempty(datacube)
        ret = 1;
    else
        ret = 0;
    end
    %% 画图
    fft_rd_map = fft(datacube,numADCSamples,1);
    fft_rd_map = abs(fft(fft_rd_map,dopplerBin_num,2));
    rd_map = sum(fft_rd_map(:,:,1:8),3);
    figure(2);imagesc(rd_map);
end