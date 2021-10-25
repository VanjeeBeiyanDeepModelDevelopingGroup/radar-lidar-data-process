%  Improved by 李嘉宝 on 13-Apr-2021
%  增加了对同一frame下不同功能段数据的区别处理

%%% This script is used to read the binary file produced by the DCA1000
%%% and Mmwave Studio
function [cpdata, vel1data, vel2data] = read5DCA1000(fid, numADCSamples, numADCSamples_vel, numChirps_cp, numChirps_vel1, numChirps_vel2, RX_num, TX_num, frame_num)
%% read file
% read .bin file

% fid = fopen(fileName,'r');
adcData = fread(fid, 'int16');
fclose(fid);
fileSize = size(adcData, 1);
% for complex data
% filesize = 2 * numADCSamples*numChirps
LVDS = zeros(1, fileSize/2);
%combine real and imaginary part into complex data
%read in file: 2I is followed by 2Q
counter = 1;
% bar = waitbar(0,'读取LVDS数据...');
for i=1:4:fileSize-1
LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
counter = counter + 2;
% if mod(i,4*10000+1) == 0 || i==fileSize-1
%     str=['读取LVDS数据...',num2str(100*i/(fileSize-1)),'%'];    % 百分比形式显示处理进程,不需要删掉这行代码就行
%     waitbar(i/(fileSize-1),bar,str);
% end
end
% close(bar);
% create column for each chirp

% cpdatanum = numADCSamples * numChirps_cp * RX_num * frame_num;
% vel1datanum = numADCSamples_vel * numChirps_vel1 * RX_num * frame_num;
% vel2datanum = numADCSamples_vel * numChirps_vel2 * RX_num * frame_num;
% LVDS_cp = reshape(LVDS(1: cpdatanum), numADCSamples, TX_num*RX_num, numChirps_cp/TX_num, frame_num);
% vel1data = permute(reshape(LVDS((cpdatanum + 1): (cpdatanum + vel1datanum)), numADCSamples_vel, RX_num * numChirps_vel1, frame_num),[2,1,3]);
% vel2data = permute(reshape(LVDS((cpdatanum + vel1datanum + 1): (cpdatanum + vel1datanum + vel2datanum)), numADCSamples_vel, RX_num * numChirps_vel2, frame_num),[2,1,3]);
% %each row is data from one chirp
% cpdata=zeros(RX_num*numChirps_cp/TX_num,numADCSamples,TX_num, frame_num);
% for i = 1:TX_num
%     cpdata(:,:,i,:)=permute(reshape(LVDS_cp(:,((i-1)*RX_num+1):((i-1)*RX_num+RX_num),:,:), numADCSamples, RX_num*numChirps_cp/TX_num,frame_num),[2,1,3]);
% end

cpdata=zeros(RX_num*numChirps_cp/TX_num,numADCSamples,TX_num, frame_num);
cpdatanum = numADCSamples * numChirps_cp * RX_num;
vel1datanum = numADCSamples_vel * numChirps_vel1 * RX_num;
vel2datanum = numADCSamples_vel * numChirps_vel2 * RX_num;
vel1data = zeros(RX_num * numChirps_vel1, numADCSamples_vel, frame_num);
vel2data = zeros(RX_num * numChirps_vel2, numADCSamples_vel, frame_num);
% bar = waitbar(0,'生成data cube...');
for frame=1:frame_num
% str=['生成data cube...',num2str(100*frame/frame_num),'%'];    % 百分比形式显示处理进程,不需要删掉这行代码就行
% waitbar(frame/frame_num,bar,str)
LVDS_slice = LVDS((1+(frame-1)*(cpdatanum+vel1datanum+vel2datanum)):(frame*(cpdatanum+vel1datanum+vel2datanum)));
vel1data(:,:,frame) = permute(reshape(LVDS_slice(1: vel1datanum), numADCSamples_vel, RX_num * numChirps_vel1),[2,1]);
vel2data(:,:,frame) = permute(reshape(LVDS_slice((vel1datanum + 1): (vel1datanum + vel2datanum)), numADCSamples_vel, RX_num * numChirps_vel2),[2,1]);
LVDS_cp = reshape(LVDS_slice((vel1datanum + vel2datanum + 1): (vel1datanum + vel2datanum + cpdatanum)), numADCSamples, TX_num*RX_num, numChirps_cp/TX_num);
% LVDS_cp = reshape(LVDS_slice(1: cpdatanum), numADCSamples, TX_num*RX_num, numChirps_cp/TX_num);
% vel1data(:,:,frame) = permute(reshape(LVDS_slice((cpdatanum + 1): (cpdatanum + vel1datanum)), numADCSamples_vel, RX_num * numChirps_vel1),[2,1]);
% vel2data(:,:,frame) = permute(reshape(LVDS_slice((cpdatanum + vel1datanum + 1): (cpdatanum + vel1datanum + vel2datanum)), numADCSamples_vel, RX_num * numChirps_vel2),[2,1]);
%each row is data from one chirp
for i = 1:TX_num
    cpdata(:,:,i,frame)=permute(reshape(LVDS_cp(:,((i-1)*RX_num+1):((i-1)*RX_num+RX_num),:,:), numADCSamples, RX_num*numChirps_cp/TX_num),[2,1]);
end 
end
% close(bar);
clear LVDS LVDS_slice LVDS_cp
end