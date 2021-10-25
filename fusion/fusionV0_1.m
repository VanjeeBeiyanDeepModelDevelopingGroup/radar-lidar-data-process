function [outConfMap,radarInLidarConfMap,new_lidarAngleGrid,new_lidarRangeGrid] = fusionV0_1(radarConfMap,lidarConfMap,radarAngleGrid,radarRangeGrid,lidarAngleGrid,lidarRangeGrid,fusionFlag)
    %% 将毫米波confMap按激光雷达confMap进行插值，得到两个大小一样的confMap
    outConfMap = [];    
    % 遍历想查询的数列，判断其范围是否在目标查询空间中，如果超出了，就去掉
    new_lidarAngleGrid = [];
    new_lidarRangeGrid = [];
    for angleIndex=1:length(lidarAngleGrid)
        angleLidar = lidarAngleGrid(angleIndex);
        if(angleLidar>=radarAngleGrid(1) && angleLidar<=radarAngleGrid(end))
            new_lidarAngleGrid = [new_lidarAngleGrid,angleLidar];
        end
    end
    for RangeIndex=1:length(lidarRangeGrid)
        rangeLidar = lidarRangeGrid(RangeIndex);
        if(rangeLidar>=lidarRangeGrid(1) && rangeLidar<=lidarRangeGrid(end))
            new_lidarRangeGrid = [new_lidarRangeGrid,rangeLidar];
        end
    end
    % 原数组
    [t,r]=meshgrid(radarAngleGrid,radarRangeGrid);
    % 想查询的数组
    [T,R]=meshgrid(new_lidarAngleGrid,new_lidarRangeGrid);
    % 插值查询，得到新的图
    radarInLidarConfMap = interp2(t,r,radarConfMap,T,R,'spline',0);
    %% 两个图像的合成
%     figure(6);
%     figure(4);imagesc(velocityCoor,distanceCoor,dopplerSum);title('单个天线的DopplerFFT波形');
%     set(gca,'YDir','normal');
%     xlabel('速度(m/s)');ylabel('距离(m)');zlabel('频谱幅值');
%     subplot(1,4,2);imagesc(new_lidarAngleGrid,new_lidarRangeGrid,radarInLidarConfMap);
%     set(gca,'YDir','normal');
%     title("毫米波雷达");
%     xlabel("角度");
%     ylabel("距离");
%     subplot(1,4,3);imagesc(lidarAngleGrid,lidarRangeGrid,lidarConfMap');
%     set(gca,'YDir','normal');
%     title("激光雷达");
%     xlabel("角度");
%     subplot(1,4,4);imagesc(lidarAngleGrid,lidarRangeGrid,outConfMap);
%     set(gca,'YDir','normal');
%     title("加和");
%     xlabel("角度");
%% 两张confMap融合，利用高斯联合概率分布计算，联合概率分布是相互独立的边缘概率分布的乘积
  % 直接取max也一样，直接取max应该是不一样的
    % outConfMap = radarInLidarConfMap+lidarConfMap';
    lidarConfMap = lidarConfMap';
    for i=1:length(new_lidarRangeGrid)
        for j=1:length(new_lidarAngleGrid)
            if fusionFlag == 0
                outConfMap(i,j) = max(radarInLidarConfMap(i,j),lidarConfMap(i,j));
            elseif fusionFlag == 1
                outConfMap(i,j) = radarInLidarConfMap(i,j)+lidarConfMap(i,j);
            elseif fusionFlag == 2
                outConfMap(i,j) = 1-(1-0.5*radarInLidarConfMap(i,j)).*(1-0.5*lidarConfMap(i,j));
            end
        end
    end
end

