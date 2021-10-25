% 目的是把所有的(u,sigma)求联合概率分布
% 高斯分布的联合概率分布是：
% u12 = u1+sigma1^2/(sigma1^2+sigma2^2)*(u2-u1)
% sigma = (sigma1^2*sigma2^2)/(sigma1^2+sigma2^2)
% 联合后生成一张置信度图用来查询存在置信度
function [confMap,allMiuSig,twoDimg,angleGrid,rangeGrid] = fusionV0_2(radarAllmiuSig,lidarAllmiuSig,radarAngleGrid,radarRangeGrid,lidarAngleGrid,lidarRangeGrid)
% 注意，输入的miusig矩阵都是“目标个数*3(angle,range,sigma)”的格式
[radarObjNum,param1] = size(radarAllmiuSig);
[lidarObjNum,param2] = size(lidarAllmiuSig);
allMiuSig = lidarAllmiuSig;
twoDimg = [];
%% 取两个雷达扫描范围的交集
minLidarRange = min(lidarRangeGrid);
maxLidarRange = max(lidarRangeGrid);
minLidarAngle = min(lidarAngleGrid);
maxLidarAngle = max(lidarAngleGrid);
minRadarRange = min(radarRangeGrid);
maxRadarRange = max(radarRangeGrid);
minRadarAngle = min(radarAngleGrid);
maxRadarAngle = max(radarAngleGrid);
% 最小测距范围
if minLidarRange <= minRadarRange
    minRange = minRadarRange;
else
    minRange = minLidarRange;
end
% 最大测距范围
if maxLidarRange <= maxRadarRange
    maxRange = maxLidarRange;
else
    maxRange = maxRadarRange;
end
% 最小测角范围
if minLidarAngle <= minRadarAngle
    minAngle = minRadarAngle;
else
    minAngle = minLidarAngle;
end
% 最大测角范围
if maxLidarAngle <= maxRadarAngle
    maxAngle = maxLidarAngle;
else
    maxAngle = maxRadarAngle;
end
% 计算离散采样空间
angleBin = 511;
rangeBin = 511;
deltaRange = (maxRange-minRange)/rangeBin;
deltaAngle = (maxAngle-minAngle)/angleBin;
rangeGrid = minRange:deltaRange:maxRange;
angleGrid = minAngle:deltaAngle:maxAngle;
confMap = zeros(rangeBin+1,angleBin+1);
if radarObjNum~=0 && lidarObjNum~=0
    %% 取接近的两个回波做联合分布
    radarPosArray = radarAllmiuSig(:,1:2);
    lidarPosArray = lidarAllmiuSig(:,1:2);
    % ns = createns(lidarPosArray,'NSMethod','kdtree');
    [idx, dist] = knnsearch(lidarPosArray,radarPosArray,'dist','euclidean');
    for i=1:length(idx)
        index_lidar = idx(i);
        if dist(i) < 5
            miu_lidar = lidarAllmiuSig(index_lidar,1:2);
            sig_lidar = lidarAllmiuSig(index_lidar,3);
            miu_radar = radarAllmiuSig(i,1:2);
            sig_radar = radarAllmiuSig(i,3);
            ksigma = sig_lidar^2/(sig_lidar^2+sig_radar^2);
            fusionMiu = miu_lidar+ksigma*(miu_radar-miu_lidar);
            fusionSig = ksigma*sig_radar^2;
%             lidarAllmiuSig(index_lidar,1:2) = fusionMiu;
            lidarAllmiuSig(index_lidar,1:2) = miu_lidar;
            lidarAllmiuSig(index_lidar,3) = fusionSig;
            allMiuSig = lidarAllmiuSig;
%         else
%             allMiuSig = [lidarAllmiuSig;radarAllmiuSig(i,:)];
        end
    end
    %% 还可以在计算激光点在毫米波的3sigma范围内，如果在，就求联合分布提升其置信度
    
end
%% 计算置信度图
[rows,cols] = size(allMiuSig);
twoDimg = zeros(size(confMap));
popIdxArray = [];
for k=1:rows
    miu = allMiuSig(k,1:2);
    sigma = allMiuSig(k,3);
    % 去掉置信度低的目标
    if sigma > 5
        popIdxArray = [k,popIdxArray];
        continue;
    end
    angle0 = floor((miu(1)-minAngle)/deltaAngle);
    range0 = floor((miu(2)-minRange)/deltaRange);
    sigma = 2*sigma;
    if angle0<=0 || range0<=0
        continue;
    end
    twoDimg(range0,angle0) = sigma;
    for m=1:rangeBin+1
        for n=1:angleBin+1
            dist = (((range0-m)*4)^2+(angle0-n)^2)/(sigma^2);
            if dist>36
                dist=36;
            end
            val = exp(-dist/2);
            if val>confMap(m,n)
                confMap(m,n)=val;
            end
        end
    end
end
% 删除置信度低的点云
for i=1:length(popIdxArray)
    idx = popIdxArray(i);
    allMiuSig(idx,:) = [];
end

