function [fusionConfMap] = fogdetector_v1(path,timeWin)
% function [fusionConfMap] = fogdetector(lidarRAmap,radarRAmap)
% 根据讨论的结果，对于正常目标点来说，回波前沿如果在远离，那么其脉宽必然是逐渐降低的
% 如果前沿逐渐靠近，那么脉宽展宽是正常的
% 但是对于雾气或者拖点来说，由于在一个光斑的视场中存在多个无法分辨的目标，会造成展宽
% 如果雾气中存在目标，回波后沿被目标决定（因为一般目标不会透光）
% 在雾气中运动的目标，对于毫米波来说，对于其测距很准确，但是角度不好确定（因为角分辨率有限）
% 但是对于激光雷达来说，一旦目标的后沿被判断出来，虽然没办法测距，因为前沿淹没在展宽里，但是测角是准确的
% 结合这两个信息，就能够从激光雷达展宽中还原出目标来了
fusionConfMap = [];
%% 存每个窗口中的5个变量，到时候要拿出来用
oldest_frontR={};
oldest_delta_frontR={};
oldest_var_frontR={};
oldest_pulseWidth={};
oldest_delta_pulseWidth={};
oldest_var_pulseWidth={};

%% 开始计算
folder = dir(path);
fileCnt = length(folder)-2;
%% 先load第一帧
filename = strcat(path,folder(3).name);
load(filename);
last_outRawFusionMap = outRawFusionMap;
%% 之后从第二帧开始
for i=4:fileCnt
    filename = strcat(path,folder(i+2).name);
    load(filename);
    %% load 到的数据
    delta_outRawFusionMap = abs(outRawFusionMap-last_outRawFusionMap);
    figure(10);
    subplot(1,2,1);imagesc(delta_outRawFusionMap);set(gca,'YDir','normal');
    subplot(1,2,2);imagesc(outRawFusionMap);set(gca,'YDir','normal');
    pause(1);
    last_outRawFusionMap = outRawFusionMap;
end
end

