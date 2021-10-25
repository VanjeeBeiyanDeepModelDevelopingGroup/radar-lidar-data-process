% 用激光来区分纯噪声与噪声干扰波形
% 对于噪声干扰波形并且有毫米波回波与之对应的目标，修正它的测距距离
% 测距角度由激光提供
function [range_angle_all] = fusionV_03(lidarADMap,radarDataCube,radarAngleGrid,radarRangeGrid,radarVelGrid,lidarAngleGrid,lidarRangeGrid)
%fusionV_03 
% 遍历每个角度下的激光回波，进行特征化
% 输出回波前沿，后沿，脉宽，评分
% 去除纯噪声波形
% 直接输出无干扰结果
% 对受干扰波形进行测距修正
% 将前沿索引去毫米波回波能量图中寻找最近的波峰
% 就是目标真实的距离
[range,angle] = size(lidarADMap);
% hungarianMat = zeros([range,angle]);
delta_radarVelGrid = radarVelGrid(2)-radarVelGrid(1);
loss = range;
range_angle_all = [];
thresh_score = 2;
thresh_pw = 30;
for i=1:angle
    range_modi_each_ang = [];
    lidarSingleLine = lidarADMap(:,i);
    [crestProp] = lidarFeatureExtract(lidarSingleLine);
    crestNum = length(crestProp)/4;
    for j=1:crestNum
        startF = crestProp(j+4*(j-1));
        endF = crestProp(j+4*(j-1)+1);
        pw = crestProp(j+4*(j-1)+2);
        score = crestProp(j+4*(j-1)+3);
        vel = crestProp(j+4*(j-1)+4);
        nearest_peakRadarIndx = 0;
        % 速度为0的时候，不要去做修正，没意义
        if vel == 0
            nearest_peakRadarIndx = startF;
        % 速度不为0时，对目标做修正
        else
            % 找到与该速度值相同的radarDataCube中的一层range-angle图
            radarRAMap_velIndx = (vel-radarVelGrid(1))/delta_radarVelGrid;
            radarRAMap = radarDataCube(:,:,radarRAMap_velIndx);
            % 取radar在目标段内进行cfar，提取一维峰值
            radarSLArr = radarRAMap(startF:endF,i);
            bit_radarSLArr = cfar_ca1D_square(radarSLArr,10,5,0.1,0);
            indx = find(bit_radarSLArr(:)==1);
            % 找到这些峰中的峰值点
            indxRadar_last = indx(1)+startF;
            radarAmp_last = radarSLArr(indxRadar_last);
            peakRadarIndx = indxRadar_last;
            maxRadarAmp = radarAmp_last;
            for k=2:length(indx)
                indxRadar_now = indx(k)+startF;
                radarAmp = radarSLArr(indxRadar_now);
                % 如果indx连续，那么说明这些值来自同一个峰，需要找到这同一个峰中的最大值
                % 好像叫：峰值合并
                if indxRadar_now-indexRadar_last == 1
                    if maxRadarAmp < radarAmp
                        peakRadarIndx = indxRadar_now;
                        maxRadarAmp = radarAmp;
                    end
                % 否则，对上一个峰值进行存储，刷新最大峰值幅值和索引    
                % 计算找到的这个毫米波峰值与激光前沿的距离，作为之后匹配的损失值，建立匈牙利权重表
                else
                    % 寻找该角度下最接近startF的毫米波峰值
                    loss_new = peakRadarIndx-startF;
                    % 更新最近峰值的index
                    if loss_new < loss
                        loss = loss_new;
                        nearest_peakRadarIndx = peakRadarIndx;
                    end
                    % 清空最大值缓存
                    peakRadarIndx = 0;
                    maxRadarAmp = 0;
                end
            end
        end
        % 通过score等值来判定要不要对这条线进行修正
        if score<thresh_score && pw>thresh_pw  % 没噪声的波，不需要被修正
            this_range = startF;
        else  % 有噪声需要被修正的波，如果没有提取到毫米波回波，该值为0，也就意味着该位置需要被修正，但是没有可对应的毫米波目标
            this_range = nearest_peakRadarIndx;
        end
        % 得到实际距离
        range_actual = lidarRangeGrid(this_range);
        angle_actual = lidarAngleGrid(i);
        % 保存出来
        range_angle_all = [range_angle_all;[angle_actual,range_actual]];
    end
end 



