function [lidarFlow] = geneLidarRangeVelMap(lidarDataFrame, radarRDMap, lidarRangeGrid,lidarAngleGrid,distanceCoor,velocityCoor)
%geneLidarRangeVelMap 计算激光的距离——速度图
% 为什么用激光的距离速度图呢，因为毫米波的距离速度图是波束收束最好，分辨率最高的
% 用光流法计算激光每个像素的速度
% 然后变换到range-velocity维度上去
%% 计算光流
% 转换激光图到笛卡尔坐标系下
[carte_lidarDataFrame,x,y] = polar2carte(lidarDataFrame, lidarAngleGrid, lidarRangeGrid, 20);
opticFlow = opticalFlowLK('NoiseThreshold',100);
lidarFlow = estimateFlow(opticFlow,carte_lidarDataFrame);
% 计算真实速度
Vx = (x(2)-x(1))/0.1*lidarFlow.Vx;
Vy = (y(2)-y(1))/0.1*lidarFlow.Vy;
figure(4);
subplot(2,2,1);imagesc(carte_lidarDataFrame);set(gca,'YDIR','normal');
hold on;
plot(lidarFlow,'DecimationFactor',[5 5],'ScaleFactor',10);
hold off
subplot(2,2,2);imagesc(x,y,Vx);set(gca,'YDIR','normal');
subplot(2,2,3);imagesc(x,y,Vy);set(gca,'YDIR','normal');hold on
% 画处速度值的大小和位置
value = min(Vy(:));
[y0,x0] = find(Vy==value);
subplot(2,2,3);text(x(x0),y(y0),num2str(value));hold on;plot(x(x0),y(y0),'rp');hold off
% subplot(2,2,4);imagesc(x,y,lidarFlow.Magnitude);set(gca,'YDIR','normal');
% 统计距离—速度二维直方图
maxVel = max(Vy(:));
minVel = min(Vy(:));
velNum = 32;
stepVel = (maxVel-minVel)/(velNum-1);
lidar_rangeDopplerMap = zeros([length(lidarRangeGrid)+1,velNum]);
for i=1:length(x)
    for j=1:length(y)
        range = sqrt(x(i)^2+y(j)^2);
        range_indx = floor(range/0.15+1);
        vel = Vy(j,i);
        vel_indx = floor((vel-minVel)/stepVel)+1;
        lidar_rangeDopplerMap(range_indx,vel_indx) = lidar_rangeDopplerMap(range_indx,vel_indx)+1; % 统计直方图
    end
end
subplot(2,2,4);imagesc(lidar_rangeDopplerMap);set(gca,'YDIR','normal');

% %% 构建该光流图对毫米波range-doppler的映射图
% % 构建lidar距离-速度图
% [lidarRangeIndx,lidarAngleIndx] = size(lidarDataFrame);
% [radarRangeIndx,radarDopplerIndex] = size(radarRDMap);
% power_lidar2radarRDMap = zeros([radarRangeIndx,radarDopplerIndex]);
% delta_radarRange = distanceCoor(2)-distanceCoor(1);
% delta_radarDoppler = velocityCoor(2)-velocityCoor(1);
% % 遍历激光range-angle图中每个点的速度（速度需要计算映射关系）
% for i = 2:lidarRangeIndx
%     for j=2:lidarAngleIndx
%         % 找到真实距离和真实速度
%         lidar_range = lidarRangeGrid(i);
% %         lidar_angle = lidarAngleGrid(j);
%         lidar_vel = lidarFlow.Vy(i,j)*0.15;
%         % 按最近邻在毫米波rnage-doppler图中寻找对应index
%         InradarRangeIndex = floor((lidar_range-distanceCoor(1))/delta_radarRange)+1;
%         InradarDopplerIndex = floor((lidar_vel-velocityCoor(1))/delta_radarDoppler)+1;
%         if lidar_range>=distanceCoor(end) || lidar_range<=distanceCoor(1) || lidar_vel >= velocityCoor(end) || lidar_vel <= velocityCoor(1)
%             continue;
%         end
%         power_lidar2radarRDMap(InradarRangeIndex,InradarDopplerIndex) = radarRDMap(InradarRangeIndex,InradarDopplerIndex);
%     end
% end
% subplot(1,2,2);imagesc(velocityCoor,distanceCoor,power_lidar2radarRDMap);set(gca,'YDIR','normal');hold off

end

