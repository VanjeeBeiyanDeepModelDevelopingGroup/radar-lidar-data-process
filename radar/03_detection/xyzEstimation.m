%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     xyzEstimation x,y,z�����Լ���Ӧ����λ�þ����ٶȼ��㺯��             %%%
%%%     maxIdx- ƽ��������������λ��                                      %%%
%%%     range- ��ѡȡ����λ�ö�Ӧ�ľ���ֵ                                   %%%
%%%     azimuthOut- ��ѡȡ����λ�õ�ˮƽ�Ƕ�fft���ֵ                       %%%
%%%     elevOut- ��ѡȡ����λ�õĴ�ֱ�Ƕ�fft���ֵ                          %%%
%%%     TX_num- ������������                                                %%%
%%%     angleBin_num- �Ƕȹ���fft bin����                                   %%%
%%%     xyzInfo- �������λ�ö�Ӧ��x,y,z����ֵ                              %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.18 version 1.1                            %%%
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