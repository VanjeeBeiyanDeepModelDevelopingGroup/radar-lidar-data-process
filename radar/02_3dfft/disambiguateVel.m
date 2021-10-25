%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     disambiguateVel function                                                                %%%
%%%     sumInputSlow- 一行Slow Chirp Sum值结果                                                  %%%
%%%     fastVel- 对应位置Fast Chirp测得的速度                                                   %%%
%%%     SumInputFast_PeakLoc- 前述步骤中Fast Chirp得到峰值位置处的Sum值                         %%%
%%%     disambVelIndx- 此处Fast Chirp峰值位置得到的速度选择（0，1 or 2）                        %%%
%%%                                                                                             %%%
%%%     Created by 李嘉宝 2021.05.25 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ peakIndx ] = disambiguateVel(sumInputSlow, fastVel, SumInputFast_PeakLoc)

    thresh = param.maxVelAssocThresh;
    MAX_VEL_ENH_NUM_NYQUIST = param.MAX_VEL_ENH_NUM_NYQUIST;
    maxUnambiguousVel = param.maxUnambiguousVel;
    invVelResolutionSlowChirp = param.invVelResolutionSlowChirp;
    validArr = zeros(1, param.N_HYPOTHESIS);
    slowChirpPeakArr = zeros(1, param.N_HYPOTHESIS);
    
    % From the fast Chirp's estimated target velocity, create a list of ambiguous velocities.
    % fastChirpAmbVel` is the ambiguous velocity for current object.
    for AmbIndx = 1: param.N_HYPOTHESIS
    
        if AmbIndx == 2
        
            threshActual = thresh;
        
        else
            % At higher velocities allow for more variation between the processed results of the 
            % 'fast chirp', and the 'slow chirp' due to range migration.
            threshActual = thresh * 2;
            
        end
        
        % Initialize the result for this ambiguous velocity to an invalid state.
        validArr(AmbIndx) = -1;
        
        % Construct the velocity hypothesis.
        AmbIndxActual = AmbIndx - 2;
        fastChirpAmbVel = (AmbIndxActual * (2 * maxUnambiguousVel)) + fastVel;

        % convert the ambiguous velocities to the slow chirp's indices.
        velIndxFlt = (fastChirpAmbVel * invVelResolutionSlowChirp);

        % Make sure that the indices lie within the slowChirp limits.
        % Also, perform a flooring operation.
        if velIndxFlt < 1
        
            velIndxTmp = int16(- velIndxFlt);

%             numMult = velIndxTmp / (2 ^ floor(log2(param.numCirpsPerFrame_vel2)));
            numMult = bitshift(velIndxTmp, -floor(log2(param.numCirpsPerFrame_vel2)), 'int16');
            velIndxFlt = (param.numCirpsPerFrame_vel2 * (numMult + 1)) + velIndxFlt;
            velIndx = floor(velIndxFlt);
        
        elseif velIndxFlt >= param.numCirpsPerFrame_vel2 + 1
            
            velIndxTmp = int16(velIndxFlt);

%             numMult = velIndxTmp / (2 ^ floor(log2(param.numCirpsPerFrame_vel2)));
            numMult = bitshift(velIndxTmp, -floor(log2(param.numCirpsPerFrame_vel2)), 'int16');
            velIndxFlt = velIndxFlt - (param.numCirpsPerFrame_vel2 * numMult);
            velIndx = floor(velIndxFlt);
        
        else
        
            velIndx = floor(velIndxFlt);
        
        end

        % Inorder to compensate for the effect of flooring, we check for both the ceil and floor of
        % slowChirpArr. i.e. we assume that either velIndx (the floor) or velIndx + 1 (the ceil),
        % could be a peak in the slowChirp's doppler array. 
        for spreadIndx = 0: param.MAX_VEL_IMPROVEMENT_NUM_SPREAD
        
            velIndxTmp = velIndx + spreadIndx;

            if velIndxTmp > param.numCirpsPerFrame_vel2
            
                velIndxTmp = velIndxTmp - param.numCirpsPerFrame_vel2;
                
            end
        
            % Is sumAbs[velIndxTmp] (or sumAbs[velIndxTmp+1]) a peak? 
            if velIndxTmp == 1
            
                prevIndx = param.numCirpsPerFrame_vel2;
            
            else
            
                prevIndx = velIndxTmp - 1;
            
            end

            if velIndxTmp == param.numCirpsPerFrame_vel2
            
                nextIndx = 1;
            
            else
            
                nextIndx = velIndxTmp + 1;
            
            end

            if (sumInputSlow(velIndxTmp) <= sumInputSlow(nextIndx)) || (sumInputSlow(velIndxTmp) <= sumInputSlow(prevIndx))
            
                continue;
                
            end
            
            if (sumInputSlow(nextIndx) == 0) || (sumInputSlow(prevIndx) == 0)
            
                continue;
            
            end
            
            % Is the amplitude of the peaks of subframe1, and subframe2 within a certain dB of each other?
            diff = abs(sumInputSlow(velIndxTmp) - SumInputFast_PeakLoc);
            
            if diff > threshActual
            
                continue;
                
            end

            % Our tests have passed, mark this as one of the possible options.
            slowChirpPeakArr(AmbIndx) = sumInputSlow(velIndxTmp);
            validArr(AmbIndx) = velIndx;

            break;
            
        end
        
    end

    % We now have a list of peak values indexed by the ambiguous velocity hypotheses.
    % Select the strongest one.
    peakIndx = -1;
    peakVal = 0;

    for AmbIndx = 1: (2 * (MAX_VEL_ENH_NUM_NYQUIST - 1)) + 1
    
        if validArr(AmbIndx) >= 0
        
            if slowChirpPeakArr(AmbIndx) > peakVal
            
                peakVal = slowChirpPeakArr(AmbIndx);
                peakIndx = AmbIndx;
            
            end
        
        end
    
    end
    
    if peakIndx ~= -1
    
        velIndx = validArr(peakIndx);
        % Since we have an association, zero out those velocity indices in the current sumAbs from
        % from being used again (in subsequent max-vel enhancement associations).
        for spreadIndx = 0: param.MAX_VEL_IMPROVEMENT_NUM_SPREAD
        
            velIndxTmp = velIndx + spreadIndx;
            if velIndxTmp > param.numCirpsPerFrame_vel2
            
                velIndxTmp = velIndxTmp - param.numCirpsPerFrame_vel2;
                
            end
            
            sumInputSlow(velIndxTmp) = 0;
            
        end
    
    else
    
        % A cheat. 
        % Always allow the detections to go through even if the disambiguation process fails. 
        % In case of failure, we use the non-disambiguated velocity.
        % The justification being that most target velocities should be below 55km/hr and 
        % that higher layer (as yet unimplemented) tracking algorithms can help disambiguate. */
        peakIndx = 2; 
     
    end

end