function [lidarDataFrames] = packLidarDataIntoFrames_v02(lidarData)
% 把雷达打包成一帧帧的形式
lidarDataFrames = [];
[rows,cols]=size(lidarData);
frame = 1;
cnt = 0;
lastLinNum = 0;
lastAngle = 0;
bar = waitbar(0,'读取激光数据...');
for i =1:rows
    if mod(i,4*1000+1) == 0 || i==rows
        str=['读取激光数据...',num2str(100*i/rows),'%'];    % 百分比形式显示处理进程,不需要删掉这行代码就行
        waitbar(i/rows,bar,str);
    end
    angle = lidarData(i,1)*256+lidarData(i,2);
%     if lidarData(i,4) == 0
%         continue
%     end
    % 11 12 13线不要，因为没几个数
    if lidarData(i,4) == 1
%         linNum = lidarData(i,3)-2+3;
        continue
    % 3 4 5线要，全是数据
    else
        linNum = lidarData(i,3)-2;
        if linNum == 1
            cnt = cnt+1;
        end
    end
    % 第4线的正反面都是以31 74开始，其他线的是以31 84开始
    linNumNow = lidarData(i,3);
%     if (angle == 31*256+74 && linNum==1) || (angle==31*256+84)
    if angle < lastAngle
        if lastLinNum ~= linNumNow
            % 帧重新开始
            frame = 0;
        end
        frame = frame+1;
        cnt = 1;
    end
    lidarDataFrames(cnt,:,linNum,frame) = lidarData(i,:); % 取到该行作为改帧第一包    
    lastLinNum = linNumNow;
    lastAngle = angle;
end
close(bar);
lidarDataFrames = lidarDataFrames(:,:,:,1:end); % 角度，距离，线序，帧号
end

