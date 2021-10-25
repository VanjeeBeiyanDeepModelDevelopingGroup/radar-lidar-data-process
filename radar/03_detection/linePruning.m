%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     linePruning function                                                                    %%%
%%%     input- 一行doppler维CFAR结果                                                            %%%
%%%     sumInput- 一行doppler sum                                                               %%%
%%%     numChirps- chirp数                                                                      %%%
%%%     MAX_NUM_DET_PER_RANGE_GATE- 每一行range gate最多峰值数量                                %%%
%%%     output- doppler维peak pruning后CFAR结果                                                 %%%
%%%     peaknum- doppler维peak pruning后CFAR结果数量                                            %%%
%%%     loc- 峰值索引位置数组                                                                   %%%
%%%                                                                                             %%%
%%%     Created by 李嘉宝 2021.05.21 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ output, peaknum, loc ] = linePruning(input, sumInput, numChirps, MAX_NUM_DET_PER_RANGE_GATE)

    peaknum = 0;
    numDetObjPerCfar = length(find(input == 1));
    cfarDetObjIndexBuf = find(input == 1);

    if numDetObjPerCfar == 0
        output = 0;
    else
        for detIdx = 1: numDetObjPerCfar

            currObjLoc = cfarDetObjIndexBuf(detIdx);

            if currObjLoc == 1

                prevIdx = numChirps;

            else

                prevIdx = currObjLoc - 1;
                
            end

            if currObjLoc == numChirps

                nextIdx = 1;

            else

                nextIdx = currObjLoc + 1;
                
            end

            if (sumInput(nextIdx) < sumInput(currObjLoc)) && (sumInput(prevIdx) < sumInput(currObjLoc))

                peaknum = peaknum + 1;

            else

                input(currObjLoc) = 0;

            end

        end

        output = input;

        % 找到最大的MAX_NUM_DET_PER_RANGE_GATE个峰值
        if MAX_NUM_DET_PER_RANGE_GATE ~= -1
            if peaknum > MAX_NUM_DET_PER_RANGE_GATE
                output = zeros(1, numChirps);
                sumInput(input == 0) = 0;
                b = sort(sumInput, 'descend');
                output(sumInput>=b(MAX_NUM_DET_PER_RANGE_GATE)) = 1;
                loc = find(output == 1);
                peaknum = length(loc);
                
            else
                loc = find(output == 1);
            end
        else
           
            loc = find(output == 1);
            
        end
        
    end
    
end