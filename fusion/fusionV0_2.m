% Ŀ���ǰ����е�(u,sigma)�����ϸ��ʷֲ�
% ��˹�ֲ������ϸ��ʷֲ��ǣ�
% u12 = u1+sigma1^2/(sigma1^2+sigma2^2)*(u2-u1)
% sigma = (sigma1^2*sigma2^2)/(sigma1^2+sigma2^2)
% ���Ϻ�����һ�����Ŷ�ͼ������ѯ�������Ŷ�
function [confMap,allMiuSig,twoDimg,angleGrid,rangeGrid] = fusionV0_2(radarAllmiuSig,lidarAllmiuSig,radarAngleGrid,radarRangeGrid,lidarAngleGrid,lidarRangeGrid)
% ע�⣬�����miusig�����ǡ�Ŀ�����*3(angle,range,sigma)���ĸ�ʽ
[radarObjNum,param1] = size(radarAllmiuSig);
[lidarObjNum,param2] = size(lidarAllmiuSig);
allMiuSig = lidarAllmiuSig;
twoDimg = [];
%% ȡ�����״�ɨ�跶Χ�Ľ���
minLidarRange = min(lidarRangeGrid);
maxLidarRange = max(lidarRangeGrid);
minLidarAngle = min(lidarAngleGrid);
maxLidarAngle = max(lidarAngleGrid);
minRadarRange = min(radarRangeGrid);
maxRadarRange = max(radarRangeGrid);
minRadarAngle = min(radarAngleGrid);
maxRadarAngle = max(radarAngleGrid);
% ��С��෶Χ
if minLidarRange <= minRadarRange
    minRange = minRadarRange;
else
    minRange = minLidarRange;
end
% ����෶Χ
if maxLidarRange <= maxRadarRange
    maxRange = maxLidarRange;
else
    maxRange = maxRadarRange;
end
% ��С��Ƿ�Χ
if minLidarAngle <= minRadarAngle
    minAngle = minRadarAngle;
else
    minAngle = minLidarAngle;
end
% ����Ƿ�Χ
if maxLidarAngle <= maxRadarAngle
    maxAngle = maxLidarAngle;
else
    maxAngle = maxRadarAngle;
end
% ������ɢ�����ռ�
angleBin = 511;
rangeBin = 511;
deltaRange = (maxRange-minRange)/rangeBin;
deltaAngle = (maxAngle-minAngle)/angleBin;
rangeGrid = minRange:deltaRange:maxRange;
angleGrid = minAngle:deltaAngle:maxAngle;
confMap = zeros(rangeBin+1,angleBin+1);
if radarObjNum~=0 && lidarObjNum~=0
    %% ȡ�ӽ��������ز������Ϸֲ�
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
    %% �������ڼ��㼤����ں��ײ���3sigma��Χ�ڣ�����ڣ��������Ϸֲ����������Ŷ�
    
end
%% �������Ŷ�ͼ
[rows,cols] = size(allMiuSig);
twoDimg = zeros(size(confMap));
popIdxArray = [];
for k=1:rows
    miu = allMiuSig(k,1:2);
    sigma = allMiuSig(k,3);
    % ȥ�����Ŷȵ͵�Ŀ��
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
% ɾ�����Ŷȵ͵ĵ���
for i=1:length(popIdxArray)
    idx = popIdxArray(i);
    allMiuSig(idx,:) = [];
end

