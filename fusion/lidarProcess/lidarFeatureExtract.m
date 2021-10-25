function [crestProp] = lidarFeatureExtract(singleLidarAD)
%lidarNoiseClassify 输出波形特征化参数，包含：
% crestPorp: 波形属性
% startF: 上升沿位置
% endF: 下降沿位置
% crestCnts: 波峰脉宽
% score: 前沿分数
% vel: 光流速度
%% 先求光流速度

%% 再求其他特征
len = length(singleLidarAD);
bit_singleLidarAD = cfar_ca1D_square_dythresh(singleLidarAD,30,15,0.2,0);
old_bit_ = bit_singleLidarAD(1);
old_val = singleLidarAD(1);
err_val = zeros(1,len);
peakStartEnd = [];
crestProp = [];
peakArray = [];
peakValue = 0;
for i=2:len
    % 求CFAR bit Map的一阶导
    new_bit_ = bit_singleLidarAD(i);
    err_bit_ = new_bit_-old_bit_;
    old_bit_ = new_bit_;
    % 求原始信号的一阶导
    new_val = singleLidarAD(i);
    err_val(i) = new_val-old_val;
    old_val = new_val;
    % 找到上升沿位置
    if err_bit_ == 1
        startF = i;
        % 求极值的位置和大小
        value = singleLidarAD(i);
        if value > peakValue
            peakValue = value;
            peakIdx = i;
        end
    % 找到下降沿位置
    elseif err_bit_ == -1
        endF = i;
        peakStartEnd = [peakStartEnd,[startF,endF]];
        startF = 0;
        endF = 0;
        % 遇到下降沿之后把极值大小和索引存起来
        peakArray = [peakArray,[peakValue,peakIdx]];
        peakValue = 0;
        peakIdx = 0;
    end
end
% 遍历波峰
crestNum = length(peakStartEnd)/2;
for crest=1:crestNum
    startF = peakStartEnd(2*crest-1);
    endF = peakStartEnd(2*crest);
    crestArray = singleLidarAD(startF:endF);
    % 求波峰脉宽
    crestCnts = length(crestArray);
    % 求波峰一阶导的最值
    err_crest = err_val(startF:endF);
    [max_err_crest,indx] = max(err_crest);
    % 计算前沿的分数
%     peakValue = peakArray(2*crest-1);
%     peakIdx = peakArray(2*crest);
    normStartF = startF/len;
    score = max_err_crest*normStartF/indx;
    if score == 0
        continue;
    end
    score = 0.5/score;
    crestProp = [crestProp;[startF,endF,crestCnts,score]];
end
        
end

