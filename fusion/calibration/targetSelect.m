function [ret] = targetSelect(rangeOut, distanceCoor, rangeBinNum, dopplerBinNum, TX_num, RX_num, lidarData_frame)
%TARGETSELECT 此处显示有关此函数的摘要
%   由于是静态标定，输入是毫米波的rangeFFT和激光的AD数据
%   从毫米波的静态rangeFFT中取到cfar结果进行channel维度的music角度分析，得到目标的角度
%   再取激光雷达AD波形中较亮的目标角度和距离
%   多次输入生成多个角度下的点对点数据，存成mat文件
%% 毫米波数据处理
%% 将输入转化成fast time* slow time* channel的形式（这是一个遗留问题）
dataCube_afterRangeFFT = zeros(rangeBinNum, dopplerBinNum, RX_num, TX_num);
for n = 1: rangeBinNum
    for i = 1: TX_num
        for m = 1: RX_num
            dataCube_afterRangeFFT(n, :, m, i) = rangeOut((0: (dopplerBinNum - 1)) * RX_num + m, n, i);
        end
    end
end
%% 做cfar，得到目标所在的距离和角度
% 取距离的平均值和角度的平均值作为这一帧的目标位置输出
L=360;
RAMap = zeros(rangeBinNum,L);
rangeOut_sum = sum(abs(rangeOut),1); % 相干累积
rangeOut_sum = log10(abs(squeeze(sum(rangeOut_sum,3))));
[biXcube] = cfar_ca1D_square(rangeOut_sum,6,4,0.1,0);
[~,cols] = find(biXcube);
finalTgtRange_radar = 0;
finalTgtAngle_radar = 0;
cnt = 0;
for i = 1:length(cols)
    % 由于角反射器固定中10m左右，故过近或者过远的index都不是我们寻找的
    if cols(i)<20 || cols(i) > 480
        continue;
    end
    AoAArr = zeros(dopplerBinNum,(TX_num-1)*RX_num);
    cnt = cnt+1;
    tgtAngle = 0;
    for j=1:dopplerBinNum
        % 在该rangebbin处取每一个chirp进行music角度分析
        AoAchannelTemp = reshape(dataCube_afterRangeFFT(cols(i),j,:,:),[1,TX_num*RX_num]);
        AoAchannelTemp = AoAchannelTemp(1:(TX_num-1)*RX_num);
%         spectAnalysis = musicAlg(AoAchannelTemp',L);
%         % 取每个chirp的角度分析峰值的平均值为角度输出结果
%         [~,index] = max(spectAnalysis);
%         tgtAngle = tgtAngle+(0.5*index-90);
        AoAArr(j,:) = AoAchannelTemp; %理论上来说，不同的chirp都应该测到相同角度的角反射器
    end
    % 直接输入该slowtime*channel数据到music算法中计算aoa
    spectAnalysis = musicAlg(AoAArr',L);
    RAMap(cols(i),:) = spectAnalysis;
    % 求峰值角度
    [~,index] = max(spectAnalysis);
    tgtAngle = (0.5*index-90);
    j = 1;
    tgtAngle = tgtAngle/j;
    finalTgtAngle_radar = finalTgtAngle_radar+tgtAngle;
    finalTgtRange_radar = finalTgtRange_radar+distanceCoor(cols(i));
end
finalTgtRange_radar = finalTgtRange_radar/cnt;
finalTgtAngle_radar = finalTgtAngle_radar/cnt;
ret = true;
if finalTgtAngle_radar > 15 || finalTgtAngle_radar < -15
    ret = false;
end
%% 激光数据处理
[pcStrc,pcPolar] = lidarRangeMeas(lidarData_frame',3);
lidarTgtAngle = pcPolar(:,1)-90;
lidarTgtRange = pcPolar(:,2);
finalTgtAngle_lidar = 0;
finalTgtRange_lidar = 0;
cnt = 0;
% 寻找照到角反上的光束，这些光束反射能量更强，AD能量更加展宽
index = find(~isnan(lidarTgtRange));
for i = 1:length(index)
    idx = index(i);
    lidarADsignal = lidarData_frame(idx,15:end);
    bi_lidarADsignal = cfar_ca1D_square(lidarADsignal,50,30,0.1,0);
%     bi_lidarADsignal = lidarADsignal > 127+20;
    gateIndx = find(bi_lidarADsignal);
    gateWidth = 0;
    allGateWidth = [];
    for j = 1:length(gateIndx)-1
        delta = gateIndx(j+1)-gateIndx(j);
        if delta == 1
            gateWidth = gateWidth+1;
        else
            % 如果门不再连续，新增该门宽度到buffer中
            allGateWidth = [allGateWidth,gateWidth];
            gateWidth = 0;
        end        
    end
    % 如果门一直连续，最后直接将该门宽度新增到buffer中
    allGateWidth = [allGateWidth,gateWidth];
    maxGateWidth = max(allGateWidth);
    if maxGateWidth > 60
        cnt = cnt +1;
        finalTgtRange_lidar = finalTgtRange_lidar+lidarTgtRange(idx);
        finalTgtAngle_lidar = finalTgtAngle_lidar+lidarTgtAngle(idx);
    end
%     figure(9);imagesc(lidarData_frame');hold on;plot(idx,0,'rp')
end
finalTgtRange_lidar = finalTgtRange_lidar/cnt;
finalTgtAngle_lidar = finalTgtAngle_lidar/cnt;
% 画图看看提取的点对不对
figure(10);
lidarData = lidarData_frame(:,40:end);
lidarData = (lidarData-min(lidarData(:)))/(max(lidarData(:))-min(lidarData(:)));
% 把两张图插值到一张图里
radarAngleGrid = (0.5*[1:L]-90);
radarRangeGrid = distanceCoor;
lidarAngleGrid = lidarTgtAngle;
lidarRangeGrid = 0.15*[1:rangeBinNum]-0.4546;
[outFusionMap,~,new_lidarAngleGrid,new_lidarRangeGrid] = fusionV0_1(RAMap,lidarData,radarAngleGrid,radarRangeGrid,lidarAngleGrid,lidarRangeGrid,1);
[lidarAngleBin,lidarRangeBin] = size(lidarData);
% figure;imagesc(lidarTgtAngle,lidarRangeGrid,lidarData');
% imagesc(new_lidarAngleGrid,new_lidarRangeGrid,outFusionMap);
imagesc(radarAngleGrid,radarRangeGrid,RAMap);title('0 deg');
set(gca,'YDIR','normal');hold on
plot(finalTgtAngle_radar,finalTgtRange_radar,'rp');hold on% xlim([-10,10]); hold on
% plot(finalTgtAngle_lidar,finalTgtRange_lidar,'bp');
result = [finalTgtAngle_radar,finalTgtRange_radar,finalTgtAngle_lidar,finalTgtRange_lidar];
% global lidar_dir_path
% index = lidar_dir_path(end-7:end-6);
% filename=['./fusion/calibration/calibData/',index,'_',datestr(now,30),'.mat'];
% save(filename,'result');
end


