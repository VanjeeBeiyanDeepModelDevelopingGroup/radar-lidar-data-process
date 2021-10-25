%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     MaxVel_2D function                                                                      %%%
%%%     inputRange- 一行Slow Chirp Sum值结果                                                    %%%
%%%     inputDoppler- 对应位置Fast Chirp测得的速度                                              %%%
%%%     sumInputFast- 前述步骤中Fast Chirp得到峰值位置处的Sum值                                 %%%
%%%     disambVelIndx- 每个Fast Chirp峰值位置计算得到的选取非模糊值（1，2，3）                  %%%
%%%     numChirps- 此处Fast Chirp峰值位置得到的速度选择（0，1 or 2）                            %%%
%%%     numDet- 测得物体数量                                                                    %%%
%%%     dopplerline- 当前处理的doppler列数                                                      %%%
%%%                                                                                             %%%
%%%     Created by 李嘉宝 2021.05.27 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ out, numDet, rangeOut_vel, velocityOut_vel, CFAROut_vel ] = MaxVel_2D(inputRange, inputDoppler, sumInputFast, disambVelIndx, ADCnum, numDet, dopplerline)

    [ out, peaknum, loc_range ] = linePruning(inputRange, sumInputFast, ADCnum, -1);
    loc_doppler = find(inputDoppler == 1);
    if peaknum > 0
        [ numDet, rangeOut_vel, velocityOut_vel, CFAROut_vel] = findandPopulateIntersectionOfDetectedObjects(peaknum, numDet, loc_range, loc_doppler, disambVelIndx, dopplerline);
    else
        numDet = numDet;
        rangeOut_vel = zeros(param.numADCSamples_vel, 1);
        velocityOut_vel = zeros(param.numADCSamples_vel, 1);
        CFAROut_vel = zeros(param.numADCSamples_vel, 1);
    end
    
end