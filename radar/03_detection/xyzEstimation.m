%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     xyzEstimation x,y,z坐标以及相应索引位置距离速度计算函数             %%%
%%%     maxIdx- 平方功率最大的索引位置                                      %%%
%%%     range- 所选取索引位置对应的距离值                                   %%%
%%%     azimuthOut- 所选取索引位置的水平角度fft结果值                       %%%
%%%     elevOut- 所选取索引位置的垂直角度fft结果值                          %%%
%%%     TX_num- 发射天线数量                                                %%%
%%%     angleBin_num- 角度估计fft bin数量                                   %%%
%%%     xyzInfo- 输出索引位置对应的x,y,z坐标值                              %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.02.18 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ xyzInfo ] = xyzEstimation( maxIdx, range, azimuthOut, elevOut, TX_num, angleBin_num)
    
    if maxIdx > (angleBin_num / 2 - 1)

        sMaxIdx = maxIdx - angleBin_num;

    else

        sMaxIdx = maxIdx;

    end

    Wx = 2 * sMaxIdx / angleBin_num;
    x = range * Wx;

    if TX_num >= 3
        peakAzimIm = imag(azimuthOut(maxIdx));
        peakAzimRe = real(azimuthOut(maxIdx));
        peakElevIm = imag(elevOut(maxIdx));
        peakElevRe = real(elevOut(maxIdx));

        Wz = atan2(peakAzimIm * peakElevRe - peakAzimRe * peakElevIm, peakAzimRe * peakElevRe + peakAzimIm * peakElevIm) / pi + (2 * Wx);

        if Wz > 1

            Wz = Wz - 2;

        elseif Wz < -1

            Wz = Wz + 2;

        end

        z = range * Wz;

    else
        
        z = 0;
        
    end
    
    y = sqrt(range * range - x * x - z * z);
    
    xyzInfo = [x, y, z];

end