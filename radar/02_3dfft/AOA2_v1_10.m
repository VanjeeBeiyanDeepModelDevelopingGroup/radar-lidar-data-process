% 生成cfar后用fft进行角度估计的RA Map结果
function [ RAMap_fft, elevOut ] = AOA2_v1_10( dopplerOut, CFAROut, TX_num, RX_num, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING)
    if isempty(CFAROut)
        azimuthOut = [];
        elevOut = [];
        return;
    end
    L = 360;
    azimuthOut = zeros(length(CFAROut(:, 1)), angleBin_num);
    elevOut = zeros(length(CFAROut(:, 1)), angleBin_num);
    % 构造RA Map矩阵
    [rangeBin_num,~,~] = size(dopplerOut);
    RAMap_fft = zeros(rangeBin_num,angleBin_num);
     for m = 1: length(CFAROut(:, 1))
         rangeIndx = CFAROut(m, 1);
         dopplerIndx = CFAROut(m, 2);
         if TX_num == 3

             Temp1 = zeros(1, angleBin_num);
             Temp2 = zeros(1, angleBin_num);

             Temp3 = reshape(dopplerOut(rangeIndx, dopplerIndx, :, :), [1, TX_num * RX_num]);
             Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
             Temp3 = dopplerCompensation(Temp3, dopplerIndx, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
             Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num));
             Temp2(1: RX_num) = Temp3((((TX_num - 1) * RX_num) + 1): (TX_num * RX_num));
             elevOut(m,:) = fftshift(fft(Temp2));
         else

             Temp1 = zeros(1, angleBin_num);
             Temp1(1: TX_num * RX_num) = reshape(dopplerOut(rangeIndx, dopplerIndx, :), [1, TX_num * RX_num]);
             Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
             Temp1 = dopplerCompensation(Temp1, dopplerIndx, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
         end
         azimuthOut(m,:) = fftshift(fft(Temp1));
         % 存进RA Map中
         RAMap_fft(rangeIndx,:) = RAMap_fft(rangeIndx,:) + fftshift(fft(Temp1));
     end
     % 这里不能取log，取log压得太狠了
%      RAMap_fft = 2*log10(abs(RAMap_fft));
    temp_Map  = abs(RAMap_fft);
    RAMap_fft = abs(RAMap_fft)/max(temp_Map(:));
end

