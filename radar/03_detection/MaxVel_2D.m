%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     MaxVel_2D function                                                                      %%%
%%%     inputRange- һ��Slow Chirp Sumֵ���                                                    %%%
%%%     inputDoppler- ��Ӧλ��Fast Chirp��õ��ٶ�                                              %%%
%%%     sumInputFast- ǰ��������Fast Chirp�õ���ֵλ�ô���Sumֵ                                 %%%
%%%     disambVelIndx- ÿ��Fast Chirp��ֵλ�ü���õ���ѡȡ��ģ��ֵ��1��2��3��                  %%%
%%%     numChirps- �˴�Fast Chirp��ֵλ�õõ����ٶ�ѡ��0��1 or 2��                            %%%
%%%     numDet- �����������                                                                    %%%
%%%     dopplerline- ��ǰ�����doppler����                                                      %%%
%%%                                                                                             %%%
%%%     Created by ��α� 2021.05.27 version 1.0                                                %%%
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