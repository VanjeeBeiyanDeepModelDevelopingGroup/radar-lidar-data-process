function [clusterFrontR,clusterPulseWidth] = lidarFnWcal(lidarRAmap,TrainWin,GuardWin,thresh)
% 计算单线扫描下每个角度的激光AD回波的前沿和脉宽
% 固定阈值，固定阈值就没办法提，因为后回波反冲后的回波，峰值都比零值低，只能cfar干
    %% 遍历激光雷达线的每个角度，计算其前沿和AD相对积分
    % 用固定阈值，好查表格
    lidarRAmap = lidarRAmap';
    [lidarLine, lidarADnum] = size(lidarRAmap);
%     thresh = 150;
    % 求每根线的前沿位置和个数
     bitMatLidarMap = [];
     for i=1:lidarLine
%          figure;plot(lidarData(i,:));
         bitMatLidarMap(i,:) = cfar_ca1D_square(lidarRAmap(i,:),TrainWin,GuardWin,thresh,0);
%          bitMatLidarMap(i,:) = cfar_ca1D_square(lidarRAmap(i,:),15,10,0.25,0);
     end
%     bilidarRAmap = lidarRAmap(:,:) > thresh; % 用固定阈值
    biLidarADmap = bitMatLidarMap;
%     biLidarADmap = CFAR2d(LidarADmap,10,10,5,5,1.2); % 用2D cfar
%     bitMatLidarMap = biLidarADmap;
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
%     [L2,num2] = bwlabel(bilidarRAmap,8); % 连通域聚类操作
%     figure;imagesc(L2');hold on
%     set(gca,'YDir','normal');
    clusterFrontR = lidarADnum*ones(lidarLine,num2); % 所有目标的前沿的index
%     rangeFrontR = zeros(lidarLine,num2); % 所有目标的前沿距离
    clusterSumAD = zeros(lidarLine,num2); % 所有目标的累计ad
    clusterSumADMax = zeros(lidarLine,num2); % 完美方波的累计ad
    clusterPulseWidth = zeros(lidarLine,num2); % 目标脉宽计数
    count = 0; % 类别计数器
    for line=1:lidarLine
        pulseWidthCnt = 0; % 脉宽计数器
        for adNum=1:lidarADnum
            % 寻找每根线的前沿
            if (L2(line,adNum)~=0)
                count = L2(line,adNum);
                pulseWidthCnt  = pulseWidthCnt+1;
                if (adNum < clusterFrontR(line,count))
                    clusterFrontR(line,count) = adNum; % 前沿
%                     rangeFrontR(line,count) = lidarRangeGrid(adNum);
                end
                clusterSumAD(line,count) = clusterSumAD(line,count)+lidarRAmap(line,adNum); % 求脉宽积分
                clusterSumADMax(line,count) = clusterSumADMax(line,count)+220; % 求完美方波脉宽积分
            else
                if pulseWidthCnt~=0
                    clusterPulseWidth(line,count) = pulseWidthCnt;
                end
                pulseWidthCnt = 0;
            end
        end
    end
    clusterSumADMax(clusterSumADMax==0)=1;
    clusterSumAD = clusterSumAD./clusterSumADMax; % 相对脉冲能量，越接近1说明波形越饱和，接近方波
%     % 画出前沿看看在哪里
%     for i=1:num2
%         plot([1:lidarLine],clusterFrontR(:,i),'go');
% %         text(clusterFrontR(:,i),[1:lidarLine],num2str(clusterSumAD(:,i),'%0.6f'),'FontSize',5);
%     end
end

