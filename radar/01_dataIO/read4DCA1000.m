%  Improved by ¿ÓºŒ±¶ on 08-Feb-2021

%%% This script is used to read the binary file produced by the DCA1000
%%% and Mmwave Studio
function [ vel1data,vel2data ] = read4DCA1000(fileName, numADCSamples_vel, numChirps_vel1, numChirps_vel2, RX_num, frame_num)

% read .bin file
fid = fopen(fileName,'r');
adcData = fread(fid, 'int16');
% if 12 or 14 bits ADC per sample compensate for sign extension
fclose(fid);
fileSize = size(adcData, 1);  %size(adcData) = 3145728 

LVDS = zeros(1, fileSize/2);
%combine real and imaginary part into complex data
%read in file: 2I is followed by 2Q
counter = 1;
for i=1:4:fileSize-1
LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
counter = counter + 2;
end

vel1datanum = numADCSamples_vel * numChirps_vel1 * RX_num;
vel2datanum = numADCSamples_vel * numChirps_vel2 * RX_num;
vel1data = zeros(RX_num * numChirps_vel1, numADCSamples_vel, frame_num);
vel2data = zeros(RX_num * numChirps_vel2, numADCSamples_vel, frame_num);

for frame=1:frame_num

LVDS_slice = LVDS((1+(frame-1)*(vel1datanum+vel2datanum)):(frame*(vel1datanum+vel2datanum)));
vel1data(:,:,frame) = permute(reshape(LVDS_slice(1: vel1datanum), numADCSamples_vel, RX_num * numChirps_vel1),[2,1]);
vel2data(:,:,frame) = permute(reshape(LVDS_slice((vel1datanum + 1): (vel1datanum + vel2datanum)), numADCSamples_vel, RX_num * numChirps_vel2),[2,1]);
    
end

end