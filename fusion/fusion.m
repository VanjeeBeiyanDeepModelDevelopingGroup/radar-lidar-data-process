function [fusionLidarADmap] = fusion(radarConfMap,lidarConfMap,radarAngleGrid,radarRangeGrid,lidarAngleGrid,lidarRangeGrid)
    %% 融合部分
    % 遍历前沿，计算位置，找到在radar图中对应的最接近目标，看是否在目标分布区域内
    % 记住：暂时是标量的融合，之后可以对此时的方差改成矩阵的
    % 毫米波的可以变成二维矩阵方差，激光雷达每根线变成一维方差（激光雷达方差似乎没什么变化）
    lidarADModifiedMat = [];
    sigmaThresh = 0.9;
    fprintf("Angle,Range,Reliability\n");
    for line=1:lidarLine
        lidarADmodifiedArray = [];
        lidarAngle = lidarAngleGrid(line); % 激光雷达目标角度
        lidarADmodifiedArray = [lidarADmodifiedArray,lidarAngle]; % 输出多线矩阵
        for frontNum=1:num2
            lidarRange = clusterFrontR(line,frontNum)/10; % 激光雷达目标距离
            lidarSigma = clusterSumAD(line,frontNum); % 激光雷达置信度
            if(lidarRange~=lidarADnum) % 该线如果没有回波，就不做操作
                % 寻找离该激光雷达前沿最近的毫米波目标
                % 先找离这个激光雷达返回值最近的毫米波图中的index
                [radarAngle,angleIndex] = min(abs(radarAngleGrid-lidarAngle));
                [radarRange,rangeIndex] = min(abs(radarRangeGrid-lidarRange));
                clusterNum = L(rangeIndex,angleIndex); % 查询该点是否在毫米波目标分布内
                % 在毫米波目标范围内的激光雷达回波
                if(clusterNum~=0) % 与毫米波位置重叠度高
                    % 找到毫米波图中这个目标的峰值位置（这里之后其实可以改成别的位置）
                    radarTargetIndex = clusterMaxPowIndex(:,clusterNum);
                    radarTargetAngle  = radarAngleGrid(radarTargetIndex(2)); % 毫米波雷达目标角度
                    radarTargetRange  = radarRangeGrid(radarTargetIndex(1)); % 毫米波雷达目标距离
                    radarSigma = clusterSumPow(clusterNum); % 毫米波雷达置信度
                    fprintf("Lidar:%3.3f, %3.3f, %3.3f\n",lidarAngle,lidarRange,lidarSigma);
                    fprintf("Radar:%3.3f, %3.3f, %3.3f\n",radarTargetAngle,radarTargetRange,radarSigma);
                    fprintf("============\n");
                    lidarADmodifiedArray = [lidarADmodifiedArray,lidarRange];  % 位置不知道怎么修改
                    lidarADmodifiedArray = [lidarADmodifiedArray,lidarSigma];  % 置信度不知道怎么修改
                elseif( abs(radarAngle-lidarAngle)<1.0 && abs(lidarRange-radarRange)<1.0) % 如果激光雷达距离与毫米波距离相差不是很远，也取可信
                    lidarADmodifiedArray = [lidarADmodifiedArray,lidarRange];  % 位置不知道怎么修改
                    lidarADmodifiedArray = [lidarADmodifiedArray,lidarSigma];  % 置信度不知道怎么修改
%                 elseif(lidarSigma>sigmaThresh) % 如果激光雷达置信度大于一定阈值，认为该回波可信
%                     lidarADmodifiedArray = [lidarADmodifiedArray,lidarRange];  % 位置不知道怎么修改
%                     lidarADmodifiedArray = [lidarADmodifiedArray,lidarSigma];  % 置信度不知道怎么修改
                else
                    lidarADmodifiedArray = [lidarADmodifiedArray,lidarRange];  % 位置直接置初始值
                    lidarADmodifiedArray = [lidarADmodifiedArray,0];  % 置信度直接置0了
                end
            end
        end
        % 计算矩阵和数列的个数，哪个短就将哪个补齐到长的长度
        lenArray = length(lidarADmodifiedArray);
        [rows,lenMatrixPerLine] = size(lidarADModifiedMat);
        if(lenMatrixPerLine<lenArray)
            addNum = abs(lenMatrixPerLine-lenArray);
            addMatrix = zeros(rows,addNum);
            lidarADModifiedMat = [lidarADModifiedMat,addMatrix];
        else
            addNum = abs(lenMatrixPerLine-lenArray);
            addMatrix = zeros(1,addNum);
            lidarADmodifiedArray = [lidarADmodifiedArray,addMatrix];
        end
        % 将结果装到矩阵中
        % 矩阵格式是：线序*[角度1，距离1，置信度1，距离2，置信度2，……]
        lidarADModifiedMat = [lidarADModifiedMat;lidarADmodifiedArray];
    end
    %% 显示结果矩阵
    [rows,cols]=size(lidarADModifiedMat);
    angleGrid = lidarADModifiedMat(:,1);
    frontMatrix = lidarADModifiedMat(:,2:end);
    for i=1:rows
        for j=1:uint8((cols-1)/2)
            reliability = frontMatrix(i,2*j);
            if(reliability~=0)
                plot(i,frontMatrix(i,2*j-1)*10,'rx');hold on
            end
        end
    end
end

