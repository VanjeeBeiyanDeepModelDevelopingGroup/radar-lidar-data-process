%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     powerAndMax �Ƕȹ������������λ��ѡ������                        %%%
%%%     azimuthOut- ˮƽ�Ƕ�AOA FFT�Ľ������Ϊ����������                   %%%
%%%     angleBin_num- angle bins����                                        %%%
%%%     maxIdx- ƽ��������������λ��                                      %%%
%%%     maxPow- �������λ�ö�Ӧ�����ƽ������ֵ                            %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.18 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ maxIdx, maxPow ] = powerAndMax( azimuthOut, angleBin_num)
    
    maxIdx = 0; 
    maxPow = 0;

    for i = 1: angleBin_num
    
        Power = power(abs(azimuthOut(i)), 2);
        if (Power > maxPow)
        
            maxPow = Power;
            maxIdx = i;
            
        end
    
    end

end