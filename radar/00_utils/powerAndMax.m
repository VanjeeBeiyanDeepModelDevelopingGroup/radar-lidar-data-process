%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     powerAndMax 角度估计最大功率索引位置选定函数                        %%%
%%%     azimuthOut- 水平角度AOA FFT的结果，作为函数的输入                   %%%
%%%     angleBin_num- angle bins数量                                        %%%
%%%     maxIdx- 平方功率最大的索引位置                                      %%%
%%%     maxPow- 输出索引位置对应的最大平方功率值                            %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.02.18 version 1.1                            %%%
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