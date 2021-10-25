function [lidarDataFrames] = packLidarDataIntoFrames(lidarData)
% 把雷达打包成一帧帧的形式
lidarDataFrames = [];
[rows,cols]=size(lidarData);
frame = 1;
cnt = 1;
for i =1:rows
    angle = lidarData(i,1)*256+lidarData(i,2);
    if angle == 23*256+132
%     if angle == 31*256+64
        % 帧开始
        frame = frame+1;
        cnt = 1;
    end
    lidarDataFrames(cnt,:,frame) = lidarData(i,:); % 取到该行作为改帧第一包    
    cnt = cnt+1;
end
lidarDataFrames = lidarDataFrames(:,:,2:end);
end

