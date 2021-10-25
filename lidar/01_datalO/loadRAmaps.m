function [lidarRAmap,radarRAmap] = loadRAmaps(frameNum,path) % loadRAmaps end
% ∂¡»°RAmap
fileCnt = num2str(frameNum);
lidar_filename = strcat(path,"lidar_",fileCnt,".mat");
radar_filename = strcat(path,"radar_",fileCnt,".mat");
lidarRAmap = load(lidar_filename);
lidarRAmap = lidarRAmap.lidarRAmap;
radarRAmap = load(radar_filename);
radarRAmap = radarRAmap.radarRAmap;
end

