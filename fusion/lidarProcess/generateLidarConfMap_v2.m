% 用新的特征生成新的激光波形置信分布图
% 其实根本的做法还是：提取->特征评价->置信评分
function [lidarConfMap,lidarMiuSig] = generateLidarConfMap_v2(LidarADmap,lidarAngleGrid,lidarRangeGrid)

[angles,ranges] = size(LidarADmap);
lidarConfMap = zeros(angles,ranges);
lidarMiuSig = [];
for j=1:angles
    singleLidarAD = LidarADmap(j,:);
    %% 取波峰
    bit_singleLidarAD = cfar_ca1D_square(singleLidarAD,50,25,0.05,0);
%     subplot(2,3,[2,3]);plot(50*bit_singleLidarAD+130);
    len = length(bit_singleLidarAD);
    old_bit_ = bit_singleLidarAD(1);
    old_val = singleLidarAD(1);
    err_val = zeros(1,len);
    peakStartEnd = [];
    startF = 1;
    endF = 1;
    peakValue = 0;
    peakIdx = 1;
    peakArray = [];
    for range=2:len
        % 求cfar的一阶导
        new_bit_ = bit_singleLidarAD(range);
        delta_bit_ = new_bit_-old_bit_;
        old_bit_ = new_bit_;
%         % 求原始信号的一阶导
%         new_val = singleLidarAD(range);
%         delta_val(range) = new_val-old_val;
%         old_val = new_val;
        % 上升沿位置
        if delta_bit_ == 1
            startF = range;
            % 求极值大小和位置
            value = singleLidarAD(range);
            if value > peakValue
                peakValue = value;
                peakIdx = range;
            end
        % 下降沿位置
        elseif delta_bit_ == -1
            endF = range;
            peakStartEnd = [peakStartEnd,[startF,endF]];
            startF = 1;
            endF = 1;
            % 遇到下降沿之后把极值大小和索引存起来
            peakArray = [peakArray,[peakValue,peakIdx]];
            peakValue = 0;
            peakIdx = 1;
        end
        value = singleLidarAD(range);
        if value > peakValue
            peakValue = value;
            peakIdx = range;
        end
    end
    crestNum = length(peakStartEnd)/2;
    % 在波峰中遍历
    for crest=1:crestNum
        startF = peakStartEnd(2*crest-1);
        endF = peakStartEnd(2*crest);
        crestArray = singleLidarAD(startF:endF);
        crestCnts = length(crestArray);
        old_val = crestArray(1);
        old_error_val = 0;
        err_err_val = zeros(1,length(singleLidarAD));
        err_val = zeros(1,length(singleLidarAD));
        peakVal = 0;
        peakIdx = 1;
        for crestCnt = 2:crestCnts
            %% 求波峰的一阶导与二阶导
            new_val = crestArray(crestCnt);
            new_error_val = new_val-old_val;
            err_val(startF+crestCnt-2) = new_error_val;
            err_err_val(startF+crestCnt-2) = new_error_val-old_error_val;
            %% 求波峰的一阶导峰值索引位置，也就是所有的二阶导零点位置（二阶导零点离散太厉害了，取不到），当然过小一阶导峰值的不要
            if new_error_val > peakVal
                peakVal = new_error_val;
                peakIdx = startF+crestCnt-2;
            end
            %% 更新
            old_val = new_val;
            old_error_val = new_error_val;
        end
        %% 计算一阶导峰值与距离门开始的位置的比作为前沿分数
%         frontScore = (peakIdx+2-startF)/crestCnts;
        startRange = lidarRangeGrid(startF)/lidarRangeGrid(end);
%         score = peakVal//(endF-peakIdx-2); % 想直接倒过来
%         score = peakVal*startRange^2/(endF-peakIdx-2); % 搞得不对呀
        score = peakVal*startRange/(peakIdx+2-startF);
        if score == 0
            continue;
        end
        sigma = 0.5/score;
%         sigma = 30*atan(score);
%         if endF-peakIdx-2 >0
        if peakIdx+2-startF > 0
            angle0 = lidarAngleGrid(j);
            range0 = lidarRangeGrid(peakIdx);
            lidarMiuSig = [angle0,range0,sigma;lidarMiuSig];% 保存angle,range,score
        %% 计算激光目标点分布图
            for m=1:angles
                for n=1:ranges
%                     sigma = sigma;
                    distance = (((j-m))^2+(peakIdx-n)^2)/(sigma^2);
                    if distance > 36 % 限定分布范围
                        distance = 36;
                    end
                    value = exp(-distance/2);
                    % 如果新算的值大于原本的值，替换掉原本的
                    if value > lidarConfMap(m,n)
                        lidarConfMap(m,n) = value;
                    end
                end
            end
        end
    end
end

end

