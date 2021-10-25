function [] = saveRAmaps(lidarRAmap,radarRAmap,frameNum,path)
% saveRAmaps 存一些.mat来用
fileCnt = num2str(frameNum);
lidar_filename = strcat(path,"lidar_",fileCnt,".mat");
radar_filename = strcat(path,"radar_",fileCnt,".mat");
save(lidar_filename,'lidarRAmap');
save(radar_filename,'radarRAmap');
% save lidar_filename lidarRAmap
% save radar_filename radarRAmap
end
