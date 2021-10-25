%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     MaxVel function                                                                         %%%
%%%     inputFast- һ��dopplerάFast Chirp����CFAR���                                          %%%
%%%     sumInputFast- һ��Fast Chirp doppler sum                                                %%%
%%%     sumInputSlow- һ��Slow Chirp doppler sum                                                %%%
%%%     numChirps- chirp��                                                                      %%%
%%%     MAX_NUM_DET_PER_RANGE_GATE- ÿһ��range gate����ֵ����                                %%%
%%%     output- dopplerάpeak pruning��CFAR���                                                 %%%
%%%     peaknum- dopplerάpeak pruning��CFAR�������                                            %%%
%%%     disambVelIndx- ÿ��Fast Chirp��ֵλ�ü���õ���ѡȡ��ģ��ֵ��1��2��3��                  %%%
%%%                                                                                             %%%
%%%     Created by ��α� 2021.05.25 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ output, peaknum, disambVelIndx ] = MaxVel(inputFast, sumInputFast, sumInputSlow, numChirps, MAX_NUM_DET_PER_RANGE_GATE)

    [ output, peaknum, loc ] = linePruning(inputFast, sumInputFast, numChirps, MAX_NUM_DET_PER_RANGE_GATE);
    CFARTHRESHOLD_N_BIT_FRAC = param.CFARTHRESHOLD_N_BIT_FRAC;
    FRAC = CFARTHRESHOLD_N_BIT_FRAC + floor(log2(param.MaxVelVirNum));
    disambVelIndx = zeros(1, param.numCirpsPerFrame_vel1);
    disambVelIndx = disambVelIndx - 2;
    
    for ele = 1: peaknum
       
        offset = quadraticInterpLog2ShortPeakLoc( sumInputFast, loc(ele), numChirps, FRAC);
        fastVelIndx = loc(ele);
%         if loc(ele) > numChirps/2
%             fastVelIndx = loc(ele) - numChirps;
%         end
        fastVelIndx = fastVelIndx - numChirps / 2;
        fastVelIndx = fastVelIndx + offset;
        fastVel = fastVelIndx * param.velResolutionFastChirp;
        [ disambVelIndx(1, loc(ele)) ] = disambiguateVel(sumInputSlow, fastVel, sumInputFast(loc(ele)));
        
    end

end