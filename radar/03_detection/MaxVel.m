%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     MaxVel function                                                                         %%%
%%%     inputFast- 一行doppler维Fast Chirp部分CFAR结果                                          %%%
%%%     sumInputFast- 一行Fast Chirp doppler sum                                                %%%
%%%     sumInputSlow- 一行Slow Chirp doppler sum                                                %%%
%%%     numChirps- chirp数                                                                      %%%
%%%     MAX_NUM_DET_PER_RANGE_GATE- 每一行range gate最多峰值数量                                %%%
%%%     output- doppler维peak pruning后CFAR结果                                                 %%%
%%%     peaknum- doppler维peak pruning后CFAR结果数量                                            %%%
%%%     disambVelIndx- 每个Fast Chirp峰值位置计算得到的选取非模糊值（1，2，3）                  %%%
%%%                                                                                             %%%
%%%     Created by 李嘉宝 2021.05.25 version 1.0                                                %%%
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