%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     peakPruning function                                                                    %%%
%%%     row- 一帧下CFAR检测出峰值的range index                                                  %%%
%%%     col- 一帧下CFAR检测出峰值的doppler index                                                %%%
%%%     peakValue- 一帧下检测出的峰值点位置幅值大小                                             %%%
%%%     numADCSamples- range bin数量                                                            %%%
%%%     dopplerBin_num- doppler bin数量                                                         %%%
%%%     peakGrpingRow- 峰值合并后的峰值结果range列表，表中为最终得出物体的 range index          %%%
%%%     peakGrpingCol- 峰值合并后的峰值结果doppler列表，表中为最终得出物体的 doppler index      %%%
%%%                                                                                             %%%
%%%     Created by 李嘉宝 2021.02.08 version 1.0                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ peakGrpingRow, peakGrpingCol ] = peakPruning(row, col, peakValue, numADCSamples, dopplerBin_num)
    peakGrpingRow = [];
    peakGrpingCol = [];
    kernel=zeros(3, 3);
    numDetectedObjects = length(row);
    peakDopplerIdx = col;
    peakRangeIdx = row;
    numDetObj = 1;

    for i = 1: numDetectedObjects

        rangeIdx = peakRangeIdx(i);
        dopplerIdx = peakDopplerIdx(i);
        kernel(2, 2) = peakValue(i);
        
%         if dopplerIdx == 16 || dopplerIdx == 17 || dopplerIdx == 18
%             continue;
%         end
        
        if rangeIdx == 1 
            kernel(1, :) = 0;
            previousRangeIdx = 0;
        else 
            previousRangeIdx = rangeIdx - 1;
        end

        if rangeIdx == numADCSamples
            kernel(3, :) = 0;
            nextRangeIdx = 0;
        else
            nextRangeIdx = rangeIdx + 1;
        end

        if dopplerIdx == 1
            leftDopplerIdx = dopplerBin_num;
        else 
            leftDopplerIdx = dopplerIdx - 1;
        end

        if dopplerIdx == dopplerBin_num
            rightDopplerIdx = 1;
        else 
            rightDopplerIdx = dopplerIdx + 1;
        end
        
        if i > 1
            if (row(i - 1) == previousRangeIdx) && (col(i - 1) == dopplerIdx)
                kernel(1, 2) = peakValue(i - 1);
            else
                kernel(1, 2) = 0;
            end
        else
            kernel(1, 2) = 0;
        end

        if i < numDetectedObjects 
            if (row(i + 1) == nextRangeIdx) && (col(i + 1) == dopplerIdx)
                kernel(3, 2) = peakValue(i + 1);
            else
                kernel(3, 2) = 0;
            end
        else
            kernel(3, 2) = 0;
        end

        if dopplerIdx == 0
            for j = numDetectedObjects: -1: 1
                if peakDopplerIdx(j) < dopplerBin_num
                    break;
                end

                if peakDopplerIdx(j) == leftDopplerIdx
                    if peakRangeIdx(j) == previousRangeIdx
                        kernel(1, 1) = peakValue(j);
                    end

                    if peakRangeIdx(j) == rangeIdx
                        kernel(2, 1) = peakValue(j);
                    end

                    if peakRangeIdx(j) == nextRangeIdx
                        kernel(3, 1) = peakValue(j);
                    end
                end
            end
        else
            for j = i: -1: 1
                if peakDopplerIdx(j) < dopplerIdx - 1 
                    break;
                end

                if peakDopplerIdx(j) == leftDopplerIdx
                    if peakRangeIdx(j) == previousRangeIdx
                        kernel(1, 1) = peakValue(j);
                    end

                    if peakRangeIdx(j) == rangeIdx
                        kernel(2, 1) = peakValue(j);
                    end

                    if peakRangeIdx(j) == nextRangeIdx
                        kernel(3, 1) = peakValue(j);
                    end
                end
            end
        end

        if dopplerIdx == dopplerBin_num
            for j = 1: numDetectedObjects
                if peakDopplerIdx(j) > 1
                    break;
                end

                if peakDopplerIdx(j) == rightDopplerIdx
                    if peakRangeIdx(j) == previousRangeIdx
                        kernel(1, 3) = peakValue(j);
                    end

                    if peakRangeIdx(j) == rangeIdx
                        kernel(2, 3) = peakValue(j);
                    end

                    if peakRangeIdx(j) == nextRangeIdx
                        kernel(3, 3) = peakValue(j);
                    end
                end
            end 
        else
            for j = i: numDetectedObjects
                if peakDopplerIdx(j) > rightDopplerIdx
                    break;
                end

                if peakDopplerIdx(j) == rightDopplerIdx
                    if peakRangeIdx(j) == previousRangeIdx
                        kernel(1, 3) = peakValue(j);
                    end

                    if peakRangeIdx(j) == rangeIdx
                        kernel(2, 3) = peakValue(j);
                    end

                    if peakRangeIdx(j) == nextRangeIdx
                        kernel(3, 3) = peakValue(j);
                    end
                end
            end
        end
        
        [kprow, kpcol] = find(kernel == max(max(kernel)));
        
        if dopplerIdx == 16 || dopplerIdx == 17 || dopplerIdx == 18
            continue;
        end
        
        if ismember([2, 2],[kprow, kpcol],'rows')
            peakGrpingRow(numDetObj) = rangeIdx;
            peakGrpingCol(numDetObj) = dopplerIdx;
            peakGrpingVal(numDetObj) = kernel(2, 2);
            numDetObj = numDetObj + 1;
        end
        
    end

end