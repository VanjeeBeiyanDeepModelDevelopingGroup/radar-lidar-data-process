function [radarConfMap,biMatRAMap,radarMiuSig] = generateRadarConfMap(radarRAmap,radarAngleGrid,radarRangeGrid)
    % 借鉴了rodnet中对于confMap的生成方式
    % 他采用了距离和类别来构造方差，之后再对confMap赋值，输出目标的位置分布
    % 公式：
    % sigma=2*arctan([class:1/2/3]/(2*rho(target)*[class_sigma:15/20/30]))
    % sigma~[class:5-15/8-20/10-30]的范围
    % 单采用毫米波生成confMap时，可以采用不一样的sigma生成方式
    % rodnet中，相当于仅仅采用了距离作为因子，那对于毫米波来说，还可以采用CFAR提取出来的波峰的相对信噪比来描述之
    %% CFAR提取信号部分
    % 毫米波雷达
    [radarRows, radarCols] = size(radarRAmap);
%     biMatRAMap = [];
%     for i=1:radarCols
%     %          figure;plot(lidarData(i,:));
%      biMatRAMap(:,i) = cfar_ca1D_square(radarRAmap(:,i),8,4,0.95,0);
%     end
    [biMatRAMap,SNRmap] = CFAR2d(radarRAmap,6,4,2,2,1.10);
%     figure;
%     imagesc(SNRmap);
%     set(gca,'YDIR','normal');
    %% 以二值图为索引，提取信号二维波形特征
    % 毫米波按照面进行处理
    % 对于毫米波来说，相当于聚类操作
    [L,num] = bwlabel(biMatRAMap,8); % 连通域聚类操作
    RGB = label2rgb(L);
    fusionLidarADmap = RGB;
%     figure;subplot(1,2,1);
%     imagesc(radarAngleGrid,radarRangeGrid,fliplr(RGB));title("所有的聚类");hold on;
%     set(gca,'YDir','normal');
    % 计算信号积分，并找到聚类的峰值代表其位置的最大估计
    clusterMaxPow = zeros(1,num); % 所有目标的最大pow
    clusterSumPow = zeros(1,num); % 所有目标的累计pow
    clusterMaxPowIndex = zeros(2,num); % 所有目标最大pow的索引
    for i=1:radarRows  % 从10开始，去掉近处的毫米波塑料盖子
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
    %% 构造confMap
    radarConfMap = zeros(radarRows, radarCols);
    radarpc = [];
    radarMiu = [radarAngleGrid(clusterMaxPowIndex(2,:));radarRangeGrid(clusterMaxPowIndex(1,:))];
    radarSig = zeros(1,num);
    for i=1:num
        rangeIndex = clusterMaxPowIndex(1,i);
        if(rangeIndex==0)
            continue
        end
        range = radarRangeGrid(rangeIndex);
        % 角度没用到，按理说也可以用到，因为毫米波电磁波并不是无差别的半空间，而是会随角度增大，能量变化
        angleIndex = clusterMaxPowIndex(2,i);
        angle = radarAngleGrid(angleIndex);
        power = clusterSumPow(i);
        % 保存目标的角度距离和能量
        radarpc(i,1) = range;
        radarpc(i,2) = angle;
        radarpc(i,3) = power;
        % 计算sigma
%         sigma = 2*atan(10/range*power)*15
        sigma = 150*atan(3*range/power)/100; % 分母越大，sigma越小，分布越集中,sigma在±30pi，我知道我为什么要把power放在分子了，因为这样符合视觉效果，如果你这个回波能量最强，
        radarSig(i) = sigma;% 那么你就分布最小，那么就看不见了，为了让你能被看见，就必须调大参数，那么其他的小能量回波就更展宽了，所以还是应该反过来，让小回波小，大回波大       
        for m=1:radarRows
            for n=1:radarCols
                distance = (((rangeIndex-m)*2)^2+(angleIndex-n)^2)/(sigma^2); % range*2 让分布在range轴变窄了
                if distance > 36 % 限定分布范围
                    distance = 36;
                end
                value = exp(-distance/2);
                % 如果新算的值大于原本的值，替换掉原本的
                if value > radarConfMap(m,n)
                    radarConfMap(m,n) = value;
                end
            end
        end
    end
    radarMiuSig = [radarMiu;radarSig]'; % 输出所有的目标的miu 和 sigma
%     figure(7);
%     subplot(1,2,1);imagesc(radarRAmap);set(gca,'YDir','normal');
%     subplot(1,2,2);imagesc(radarConfMap);title("所有的聚类");hold on;
%     set(gca,'YDir','normal');
%         % 画出目标的位置和积分值
%     subplot(1,2,2);plot(clusterMaxPowIndex(2,:),clusterMaxPowIndex(1,:),'ro');
%     subplot(1,2,2);text(clusterMaxPowIndex(2,:),clusterMaxPowIndex(1,:),num2str((clusterSumPow(:)),'%0.6f'),'FontSize',10);
end