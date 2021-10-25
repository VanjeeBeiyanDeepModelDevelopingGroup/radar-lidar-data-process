%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     quadraticInterpLog2ShortPeakLoc function                            %%%
%%%     locsbuffer- 当前range bin位置的一行doppler数据                      %%%
%%%     indx- 之前步骤得到的峰值位置索引                                    %%%
%%%     len- doppler bin数量                                                %%%
%%%     frac- 所占小数位                                                    %%%
%%%     offset- 计算速度可以做的精细化索引偏移                              %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.05.25 version 1.0                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ offset ] = quadraticInterpLog2ShortPeakLoc( locsbuffer, indx, len, frac)

    y0 = antilog2(locsbuffer(indx), frac); 
    
    % Circular shifting for finding the neighbour at the extremes.
    if indx == len
        yp1 = antilog2(locsbuffer(1), frac); 
    else
        yp1 = antilog2(locsbuffer(indx + 1), frac);
    end
    
    if indx == 1
        ym1 = antilog2(locsbuffer(len), frac); 
    else
        ym1 = antilog2(locsbuffer(indx - 1), frac);
    end
    
    Den = (2 * ((2 * y0) - yp1 - ym1));
    
    % A reasonable restriction on the inverse. 
    % Note that y0 is expected to be larger than 
    % yp1 and ym1. 
    if (Den > 0.15)
        % Compute the interpolated maxima as 
        %   thetaPk = (yp1 - ym1)/(2 * ((2 * y0) - yp1 - ym1))
        offset = (yp1 - ym1)/Den;
    else
        offset = 0;
    end

end