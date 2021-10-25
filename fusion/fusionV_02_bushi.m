%这个函数的目的是，直接弄一张画布，把confmap画在这个画布上
% 借鉴了rodnet中对于confMap的生成方式
% 他采用了距离和类别来构造方差，之后再对confMap赋值，输出目标的位置分布
% 公式：
% sigma=2*arctan([class:1/2/3]/(2*rho(target)*[class_sigma:15/20/30]))
% sigma~[class:5-15/8-20/10-30]的范围
% 单采用毫米波生成confMap时，可以采用不一样的sigma生成方式
% rodnet中，相当于仅仅采用了距离作为因子，那对于毫米波来说，还可以采用CFAR提取出来的波峰的相对信噪比来描述之

%%%%%%%%%%%%%%%%%%%%%%%%%没有用%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
function [outConfMap] = fusionV_02(LidarADmap,lidarAngleGrid,lidarRangeGrid,radarRAmap,radarAngleGrid,radarRangeGrid)
    %% 寻找radar中的目标，计算均值和方差
    %% CFAR提取信号部分
    % 毫米波雷达
    [radarRows, radarCols] = size(radarRAmap);
    biMatRAMap = CFAR2d(radarRAmap, 1.15);
    %% 以二值图为索引，提取信号二维波形特征
    % 毫米波按照面进行处理
    % 对于毫米波来说，相当于聚类操作
    [L,num] = bwlabel(biMatRAMap,8); % 连通域聚类操作
    RGB = label2rgb(L);
    fusionLidarADmap = RGB;
    % 计算信号积分，并找到聚类的峰值代表其位置的最大估计
    clusterMaxPow = zeros(1,num); % 所有目标的最大pow
    clusterSumPow = zeros(1,num); % 所有目标的累计pow
    clusterMaxPowIndex = zeros(2,num); % 所有目标最大pow的索引
    for i=10:radarRows  % 从10开始，去掉近处的毫米波塑料盖子
        for j=1:radarCols
            if(L(i,j)~=0)
                count = L(i,j);
                pow = radarRAmap(i,j);
                if pow>clusterMaxPow(count)
                    clusterMaxPow(count) = pow;
                    clusterMaxPowIndex(1,count) = i;
                    clusterMaxPowIndex(2,count) = j;
                end
                clusterSumPow(count) = clusterSumPow(count)+pow;
            end
        end
    end
%     clusterSumPow = clusterSumPow/max(clusterSumPow);
    clusterSumPow = log10(clusterSumPow);
    %% 寻找lidar中的目标，计算均值和方差
    % 遍历激光雷达线的每个角度，计算其前沿和AD相对积分
    % 用固定阈值，好查表格
    [lidarLine, lidarADnum] = size(LidarADmap);
    thresh = 150;
    % 求每根线的前沿位置和个数
     bitMatLidarMap = [];
     for i=1:lidarLine
%          bitMatLidarMap(i,:) = cfar_ca1D_square(LidarADmap(i,:),15,10,0.02,0);
         bitMatLidarMap(i,:) = cfar_ca1D_square(LidarADmap(i,:),15,10,0.3,0);
     end
%     biLidarADmap = LidarADmap(:,:) > thresh; % 用固定阈值
    biLidarADmap = bitMatLidarMap;
%     biLidarADmap = bitMatLidarMap;
    % 对每条线进行聚类，对波峰计数编号
    L2 = zeros(lidarLine, lidarADnum);
    for j=1:lidarLine
        num2 = 0;
        neighborIndex = 0;
        for i=1:lidarADnum
            if(biLidarADmap(j,i)~=0)
                if(abs(i-neighborIndex)==1)
                    L2(j,i) = num2;
                else
                    num2 = num2+1;
                    L2(j,i) = num2;
                end
                neighborIndex = i;
            end
        end
    end
    num2 = max(L2(:));
    clusterFrontR = lidarADnum*ones(lidarLine,num2); % 所有目标的前沿ad
    clusterSumAD = zeros(lidarLine,num2); % 所有目标的累计ad
    clusterSumADMax = zeros(lidarLine,num2); % 完美方波的累计ad
    for line=1:lidarLine
        for adNum=1:lidarADnum
            % 寻找每根线的前沿
            if (L2(line,adNum)~=0)
                count = L2(line,adNum);
                if (adNum < clusterFrontR(line,count))
                    clusterFrontR(line,count) = adNum; % 前沿
                end
                clusterSumAD(line,count) = clusterSumAD(line,count)+LidarADmap(line,adNum); % 求脉宽积分
                clusterSumADMax(line,count) = clusterSumADMax(line,count)+220; % 求脉宽积分
            end
        end
    end
    clusterSumADMax(clusterSumADMax==0)=1;
    clusterSumAD = clusterSumAD./clusterSumADMax; % 相对积分比，越接近1说明波形越饱和，接近方波

    
    %% 大家都用lidar的画布，因为lidar有更小的视场，更长的距离
    % 那么这里其实就是需要把radar的(r,theta)映射到lidar坐标系(r,theta)下去，就可以了
    % 应有的操作是：1.空间位置标定 2.插值操作
    % 在这里由于只有azimuth的，所以平移只有左右没有上下，但是有旋转
%     clusterFrontR = lidarADnum*ones(lidarLine,num2); % 所有目标的前沿ad
%     clusterSumAD = zeros(lidarLine,num2); % 所有目标的累计ad
    radarInLidar_R = [];
    radarInLidar_Theta = [];
    radarR_array = [];
    radarTheta_array = [];
    %% 这里是标定两个坐标的变换和平移
    for i=1:num
        radarRindex = clusterMaxPowIndex(1,i);
        radarThetaIndex = clusterMaxPowIndex(2,i);
        if radarRindex<=0 || radarThetaIndex<=0
            continue
        end
        radarR = radarRangeGrid(radarRindex);
        radarTheta = radarAngleGrid(radarThetaIndex)*pi/180;
        radarR_array = [radarR_array,radarR];
        radarTheta_array = [radarTheta_array,radarTheta*180/pi];
        % 以中间90度为起始0点
        radarX = radarR*sin(radarTheta);
        radarY = radarR*cos(radarTheta);
        radarZ = 1;
        % 旋转平移到lidar坐标系
        % 假设旋转矩阵为R=diag(3) 平移矩阵为t=[0,0,0];
        R = diag([1,1,1],0);
        t = [0;0;0];
        radarINlidarPos = R*[radarX;radarY;radarZ]+t;
        radarXinlidar = radarINlidarPos(1);
        radarYinlidar = radarINlidarPos(2);
        radarZinlidar = radarINlidarPos(3);
        % 转回极坐标
        % 以中间90度为起始0点
        theta = atan(radarXinlidar/radarYinlidar)*180/pi;
        r = sqrt(radarXinlidar^2+radarYinlidar^2);
        radarInLidar_R = [radarInLidar_R,r];
        radarInLidar_Theta = [radarInLidar_Theta,theta];
    end
%     [t,r] = meshgrid(radarAngleGrid,radarRangeGrid);
%     [T,R] = meshgrid(radarInLidar_X,radarInLidar_Y);
    % 看看变换出来的点在不在附近
    figure;
    subplot(1,2,1);imagesc(radarAngleGrid,radarRangeGrid,radarRAmap);
    set(gca,'YDir','normal');hold on
    plot(radarInLidar_Theta,radarInLidar_R,'ro');hold on
    plot(radarTheta_array,radarR_array,'gx');
%     text(radarTheta_array,radarR_array,num2str((clusterMaxPow),'%0.6f'),'FontSize',10);
    subplot(1,2,2);imagesc(lidarAngleGrid,lidarRangeGrid,LidarADmap');
    set(gca,'YDir','normal');hold on
    plot(radarInLidar_Theta,radarInLidar_R,'ro');
    
    %% 现在把radar目标位置对应的lidar索引位置插值出来
    % 最近邻查找吧
    
    %% 现在在lidar画布上生成confMap
    
end

