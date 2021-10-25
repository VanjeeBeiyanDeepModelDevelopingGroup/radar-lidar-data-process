function [pcdPerFrame] = outputpointcloud(fusionConfMap,confAngleGrid,confRangeGrid,lidar_pc,radar_pc)
% 本函数输出每帧所有的点云
% radarpc: range,angle,power - 2d -插值之前的位置
% lidarpc: frontR,pulseWidth,sumADratio - 3d
% 查找输入点位置的置信度
% 输出格式为：x,y,z,intensity(脉宽),reliability
pcdPerFrame = [];
%% 构建坐标轴kdtree以供查找
% tree = createns(confRangeGrid,'NSMethod','kdtree');
%% 激光雷达点云置信度查找
lidarPCfrontR = lidar_pc(:,:,1);
[angles,numR] = size(lidarPCfrontR);
reliab = [];
figure(5);
imagesc(fusionConfMap);
set(gca,'YDir','normal');hold on;
for i=1:angles
    for j=1:numR
        frontR = lidarPCfrontR(i,j);
        step = confRangeGrid(2)-confRangeGrid(1);
        idx = floor((frontR-confRangeGrid(1))/step)+1;
        if idx <= 0
            reliab = [reliab,0];
            continue
        end
        reliab = [reliab,fusionConfMap(idx,i)];
%         plot(idx,i,'rx');
        plot(i,idx,'rx');hold on
%         text(i,idx,num2str(fusionConfMap(idx,i))); hold on
    end
end
lidar_pc(:,:,4) = reliab;
end

