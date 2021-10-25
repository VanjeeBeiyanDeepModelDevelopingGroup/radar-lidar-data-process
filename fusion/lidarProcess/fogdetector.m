function [fusionConfMap] = fogdetector(path,timeWin)
% function [fusionConfMap] = fogdetector(lidarRAmap,radarRAmap)
% 根据讨论的结果，对于正常目标点来说，回波前沿如果在远离，那么其脉宽必然是逐渐降低的
% 如果前沿逐渐靠近，那么脉宽展宽是正常的
% 但是对于雾气或者拖点来说，由于在一个光斑的视场中存在多个无法分辨的目标，会造成展宽
% 如果雾气中存在目标，回波后沿被目标决定（因为一般目标不会透光）
% 在雾气中运动的目标，对于毫米波来说，对于其测距很准确，但是角度不好确定（因为角分辨率有限）
% 但是对于激光雷达来说，一旦目标的后沿被判断出来，虽然没办法测距，因为前沿淹没在展宽里，但是测角是准确的
% 结合这两个信息，就能够从激光雷达展宽中还原出目标来了
fusionConfMap = [];
%% 维度统一参数
last_angles_rows = 0;
last_frontNum_cols = 0;
now_angles_rows = 0;
now_frontNum_cols = 0;
%% 存每个窗口中的5个变量，到时候要拿出来用
oldest_frontR={};
oldest_delta_frontR={};
oldest_var_frontR={};
oldest_pulseWidth={};
oldest_delta_pulseWidth={};
oldest_var_pulseWidth={};
%% 脉冲前沿的变量
last_frontR = [];
new_frontR = [];
sum_frontR = [];
sum_delta_frontR = [];
sum_var_frontR = [];
ave_frontR = [];
ave_delta_frontR = [];
var_frontR = [];
%% 脉宽的变量
last_pulseWidth = [];
new_pulseWidth = [];
sum_pulseWidth = [];
sum_delta_pulseWidth = [];
sum_var_pulseWidth = [];
ave_pulseWidth = [];
ave_delta_pulseWidth = [];
var_pulseWidth = [];
%% 计数器
timeWinCnt = 0;
outputCnt = 0;
for i=1:100 % 从第二帧开始
    %% load数据
    [lidarRAmap,radarRAmap] = loadRAmaps(i,path);
    [rangeIdxFrontR,clusterPulseWidth] = lidarFnWcal(lidarRAmap,15,10,0.25);
    %% 把相邻帧的数据维度统一，还有累加器的维度也要统一，是我写的太麻烦还是matlab就是麻烦
    [now_angles_rows,now_frontNum_cols] = size(rangeIdxFrontR);
    % 当前帧比上一帧维度大
    if now_frontNum_cols>last_frontNum_cols
        % 更新维度之后要更新维度大小（如果大是不需要的，因为都是小的补到大的去）
%         now_frontNum_cols = last_frontNum_cols;
        % 前沿
        last_frontR_fill = 564*ones(now_angles_rows,now_frontNum_cols);
        last_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = last_frontR;
        last_frontR = last_frontR_fill;
        % 脉宽
        last_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        last_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = last_pulseWidth;
        last_pulseWidth = last_pulseWidth_fill;
        % 前沿位置累加器：sum_frontR
        sum_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_frontR;
        sum_frontR = sum_frontR_fill;
        % 前沿变化率累加器：sum_delta_frontR
        sum_delta_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_delta_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_delta_frontR;
        sum_delta_frontR = sum_delta_frontR_fill;
        % 前沿方差累加器：sum_var_frontR
        sum_var_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_var_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_var_frontR;
        sum_var_frontR = sum_var_frontR_fill;
        % 脉宽大小累加器：sum_pulseWidth
        sum_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_pulseWidth;
        sum_pulseWidth = sum_pulseWidth_fill;
        % 脉宽变化率累加器：sum_delta_pulseWidth
        sum_delta_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_delta_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_delta_pulseWidth;
        sum_delta_pulseWidth = sum_delta_pulseWidth_fill;
        % 脉宽方差累加器：sum_var_pulseWidth
        sum_var_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_var_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_var_pulseWidth;
        sum_var_pulseWidth = sum_var_pulseWidth_fill;
        % 前沿位置均值
        ave_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        ave_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = ave_frontR;
        ave_frontR = ave_frontR_fill;
        % 脉宽大小均值
        ave_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        ave_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = ave_pulseWidth;
        ave_pulseWidth = ave_pulseWidth_fill;
    % 当前帧比上一帧维度小
    elseif now_frontNum_cols<last_frontNum_cols
        % 前沿
        rangeIdxFrontR_fill = 564*ones(last_angles_rows,last_frontNum_cols);
        rangeIdxFrontR_fill(1:now_angles_rows,1:now_frontNum_cols) = rangeIdxFrontR;
        rangeIdxFrontR = rangeIdxFrontR_fill;
        % 脉宽
        clusterPulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
        clusterPulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = clusterPulseWidth;
        clusterPulseWidth = clusterPulseWidth_fill;
%         % 前沿位置累加器：sum_frontR
%         sum_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_frontR;
%         sum_frontR = sum_frontR_fill;
%         % 前沿变化率累加器：sum_delta_frontR
%         sum_delta_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_delta_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_delta_frontR;
%         sum_delta_frontR = sum_delta_frontR_fill;
%         % 前沿方差累加器：sum_var_frontR
%         sum_var_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_var_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_var_frontR;
%         sum_var_frontR = sum_var_frontR_fill;
%         % 脉宽大小累加器：sum_pulseWidth
%         sum_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_pulseWidth;
%         sum_pulseWidth = sum_pulseWidth_fill;
%         % 脉宽变化率累加器：sum_delta_pulseWidth
%         sum_delta_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_delta_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_delta_pulseWidth;
%         sum_delta_pulseWidth = sum_delta_pulseWidth_fill;
%         % 脉宽方差累加器：sum_var_pulseWidth
%         sum_var_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_var_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_var_pulseWidth;
%         sum_var_pulseWidth = sum_var_pulseWidth_fill;
%         % 前沿位置均值
%         ave_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         ave_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = ave_frontR;
%         ave_frontR = ave_frontR_fill;
%         % 脉宽大小均值
%         ave_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         ave_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = ave_pulseWidth;
%         ave_pulseWidth = ave_pulseWidth_fill;
        % 更新维度之后要更新维度大小（如果大是不需要的，因为都是小的补到大的去）
        now_frontNum_cols = last_frontNum_cols;
    % 如果相同就直接跳过了
    end
    %% 前沿
    new_frontR = rangeIdxFrontR;
    sum_frontR = sum_frontR+new_frontR;
    sum_delta_frontR = sum_delta_frontR+(new_frontR-last_frontR);
    sum_var_frontR = sum_var_frontR+abs(new_frontR-ave_frontR);
    % 缓存
    oldest_frontR{end+1} = new_frontR;
    oldest_delta_frontR{end+1} = new_frontR-last_frontR;
    oldest_var_frontR{end+1} = abs(new_frontR-ave_frontR);
    %% 脉宽
    new_pulseWidth = clusterPulseWidth;
    sum_pulseWidth = sum_pulseWidth+new_pulseWidth;
    sum_delta_pulseWidth = sum_delta_pulseWidth+(new_pulseWidth-last_pulseWidth);
    sum_var_pulseWidth = sum_var_pulseWidth+abs(new_pulseWidth-ave_pulseWidth);
    % 缓存
    oldest_pulseWidth{end+1} = new_pulseWidth;
    oldest_delta_pulseWidth{end+1} = new_pulseWidth-last_pulseWidth;
    oldest_var_pulseWidth{end+1} = abs(new_pulseWidth-ave_pulseWidth);
    %% 计数器（没意义）
%     timeWinCnt = timeWinCnt+1;
    if i > timeWin
        %% 计算关键特征参数
        ave_frontR = sum_frontR/timeWin; % 前沿均值位置
        ave_delta_frontR = sum_delta_frontR/timeWin; % 前沿窗口均值变化率
        var_frontR = sum_var_frontR/timeWin; % 前沿方差
        
        ave_pulseWidth = sum_pulseWidth/timeWin; % 脉宽均值
        ave_delta_pulseWidth = sum_delta_pulseWidth/timeWin;% 脉宽变化率
        var_pulseWidth = sum_var_pulseWidth/timeWin; % 脉宽方差
        
        %% 限制累加变量
        sum_frontR  = sum_frontR-oldest_frontR{1};
        sum_delta_frontR = sum_delta_frontR-oldest_delta_frontR{1};
        sum_var_frontR = sum_var_frontR-oldest_var_frontR{1};
        
        sum_pulseWidth = sum_pulseWidth-oldest_pulseWidth{1};
        sum_delta_pulseWidth = sum_delta_pulseWidth-oldest_delta_pulseWidth{1};
        sum_var_pulseWidth = sum_var_pulseWidth-oldest_var_pulseWidth{1};
        
        %% 清空掉最老的那个
        oldest_frontR(1)=[];
        oldest_delta_frontR(1)=[];
        oldest_var_frontR(1)=[];
        oldest_pulseWidth(1)=[];
        oldest_delta_pulseWidth(1)=[];
        oldest_var_pulseWidth(1)=[];
    end
    %% 更新
    last_frontR = new_frontR;
    last_pulseWidth = new_pulseWidth;
    last_angles_rows=now_angles_rows;
    last_frontNum_cols=now_frontNum_cols;
    %% 进行特征判断
    % ave_frontR = sum_frontR/timeWinCnt; % 前沿均值位置
    % ave_delta_frontR = sum_delta_frontR/timeWinCnt; % 前沿窗口均值变化率
    % var_frontR = sum_var_frontR/timeWinCnt; % 前沿方差        
    % ave_pulseWidth = sum_pulseWidth/timeWinCnt; % 脉宽均值
    % ave_delta_pulseWidth = sum_delta_pulseWidth/timeWinCnt;% 脉宽变化率
    % var_pulseWidth = sum_var_pulseWidth/timeWinCnt; % 脉宽方差
    if isempty(ave_frontR) || isempty(ave_delta_frontR) || isempty(var_frontR) || isempty(ave_pulseWidth) || isempty(ave_delta_pulseWidth) || isempty(var_pulseWidth)
        fprintf("还才第%d帧，现在还没有输出\n",i);
        continue
    end
    [angleRow,peakNum] = size(ave_frontR); % 用哪个都一样，维度都是一样的
%     for angle=1:angleRow
%         for num=1:peakNum
%             %% 前沿
    
    
end
end

