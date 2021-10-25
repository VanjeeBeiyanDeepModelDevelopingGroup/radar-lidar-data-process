%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     findandPopulateIntersectionOfDetectedObjects function                                   %%%
%%%     sumInputSlow- һ��Slow Chirp Sumֵ���                                                  %%%
%%%     fastVel- ��Ӧλ��Fast Chirp��õ��ٶ�                                                   %%%
%%%     SumInputFast_PeakLoc- ǰ��������Fast Chirp�õ���ֵλ�ô���Sumֵ                         %%%
%%%     disambVelIndx- �˴�Fast Chirp��ֵλ�õõ����ٶ�ѡ��0��1 or 2��                        %%%
%%%                                                                                             %%%
%%%     Created by ��α� 2021.05.27 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ numDetObj2D, rangeOut_vel, velocityOut_vel, CFAROut_vel ] = findandPopulateIntersectionOfDetectedObjects(numDetObjPerCfar, numDetObj2D, loc, loc_doppler, disambVelIndx, dopplerline)

%     oneQFormat = bitshift(1, param.xyzOutputQFormat);
    rangeOut_vel = zeros(param.numADCSamples_vel, 1);
    velocityOut_vel = zeros(param.numADCSamples_vel, 1);
    CFAROut_vel = zeros(param.numADCSamples_vel, 1);
    
    % For each object in the CFAR detected object list
    for detIdx2 = 1: numDetObjPerCfar
    
        % if there is space in the detObj2DRaw matrix
        if numDetObj2D < param.maxNumObj2DRaw
        
            % locate the 1D CFAR corresponding to the current range gate
            
            rangeIdx = loc(detIdx2);

            % Check if the 1D CFAR matches the current objects doppler. 
            % Also, check if the velocity disambiguation output is valid.
            if (disambVelIndx(rangeIdx) >= 0) && ismember(rangeIdx, loc_doppler)

                % Calculate 
                % 1. The speed (after disambiguation). 
                if dopplerline > (param.numCirpsPerFrame_vel1 / 2)

                    dopplerIdxActual = dopplerline - param.numCirpsPerFrame_vel1;

                else

                    dopplerIdxActual = dopplerline;

                end
               
                disambiguatedSpeed = (dopplerIdxActual * param.velResolutionFastChirp) + (2 * (disambVelIndx(rangeIdx) - 2)) * param.maxUnambiguousVel;

                % 2. The range (after correcting for the antenna delay).
                %    MIN_RANGE_OFFSET_METERS- The radar's range estimate has a constant error due to the finite distance from the antenna to the LO.
                range = (rangeIdx * param.rangeResolution) - param.MIN_RANGE_OFFSET_METERS;

                if (range < 0)

                    range = 0;

                end

                % 3. Populate the output Array.
                if disambVelIndx(rangeIdx) == 2
                CFAROut_vel(rangeIdx) = 1;
                velocityOut_vel(rangeIdx) =  disambiguatedSpeed; 
                rangeOut_vel(rangeIdx) = range;
                
                numDetObj2D = numDetObj2D + 1;
                end
            end
        end
    end
end
         