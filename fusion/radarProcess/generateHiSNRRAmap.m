function [ramap] = generateHiSNRRAmap(dopplerOut,azimuthout)
% 生成高信噪比range-angle heat map
% 原理是：视场中有速度的目标，他的range-angle heat map信噪比更高，这样的话，取在目标速度附近的
% RA-heat map叠加（或者取最大），就可以得到更高信噪比的ra-map

% 在速度维上取1D-CFAR
    %% 获取rangedopplermap
    dopplerLog2Abs = log2(abs(dopplerOut));
    dopplerSum = sum(dopplerLog2Abs, [3 4]);
    rangedopplermap = squeeze(dopplerSum); 
    %% CFAR提取doppler图
    [radarRows, radarCols] = size(rangedopplermap);
    biMatRDMap = [];
    for i=1:radarRows
    %          figure;plot(lidarData(i,:));
     biMatRDMap(i,:) = cfar_ca1D_square(rangedopplermap(i,:),8,4,0.25,0);
    end
    %% 以二值图为索引，提取信号二维波形特征
    % 毫米波按照面进行处理
    % 对于毫米波来说，相当于聚类操作
    [L,num] = bwlabel(biMatRDMap,8); % 连通域聚类
%     RGB = label2rgb(L);
%     figure;
%     imagesc(RGB);title("所有的聚类");hold on;
%     set(gca,'YDir','normal');
    colScale = [zeros(num,1),radarCols*ones(num,1)];
    rowScale = zeros(num,2);
    for i=1:radarRows
        for j=1:radarCols
            if L(i,j)~=0
                classNum = L(i,j);
                colNow = j;
                % 找最大
                if colNow > colScale(classNum,1)
                    colScale(classNum,1) = colNow;
                    rowScale(classNum,1) = i;
                end
                % 找最小
                if colNow < colScale(classNum,2)
                    colScale(classNum,2) = colNow;
                    rowScale(classNum,2) = i;
                end
            end
        end
    end
%             
%     figure(10);
%     imagesc(biMatRDMap);hold on
%     set(gca,'YDir','normal');
%         % 画出目标的位置和积分值
%     plot(colScale(:,1),rowScale(:,1),'ro');hold on
%     plot(colScale(:,2),rowScale(:,2),'bo');
% %     text(clusterMaxPowIndex(2,:),clusterMaxPowIndex(1,:),num2str((clusterSumPow(:)),'%0.6f'),'FontSize',10);
    %% 在azimuthout中挑选colScale中有的速度页，加起来
     % 取操作
    [ranges,antennas,dopplers] = size(azimuthout); % 64个antenna采样值，因为是补0的
    indexArray = zeros(1,dopplers);
    newAzimuthSum = ones(ranges,antennas);    
    for i=1:num
        startCol = colScale(i,2);
        endCol = colScale(i,1);
        for j=startCol:endCol
            indexArray(j) = 1;
        end
    end
%     indexArray
%     figure(11);
    cnt = 0;
    for dop=1:dopplers
        if indexArray(dop) == 1
            if abs(dop-dopplers/2) < 3
                continue
            end
%             fprintf("取的帧号：%d\n",dop);
            cnt = cnt+1;
%             newAzimuthSum = newAzimuthSum+azimuthout(:,:,dop);
            newAzimuthSum = newAzimuthSum.*azimuthout(:,:,dop);
%             subplot(4,8,dop);imagesc(azimuthout(:,:,dop));
%             set(gca,'YDir','normal');axis off
        end
    end
%     indexArray
    ramap = newAzimuthSum;
%     figure(12);
%     imagesc(newAzimuthSum);hold on
%     set(gca,'YDir','normal');
end

