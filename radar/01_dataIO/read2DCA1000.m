%  Improved by Àî¼Î±¦ on 08-Feb-2021

%%% This script is used to read the binary file produced by the DCA1000
%%% and Mmwave Studio
function [retVal] = read2DCA1000(fileName, numADCSamples, RX_num, TX_num, frame_num)
%% global variables
% change based on sensor config
numADCBits = 16; % number of ADC bits per sample
isReal = 0; % set to 1 if real only data, 0 if complex data0
%% read file
% read .bin file
fid = fopen(fileName,'r');
adcData = fread(fid, 'int16');
% if 12 or 14 bits ADC per sample compensate for sign extension
if numADCBits ~= 16
l_max = 2^(numADCBits-1)-1;
adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
end
fclose(fid);
fileSize = size(adcData, 1);  %size(adcData) = 3145728 
% real data reshape, filesize = numADCSamples*numChirps
if isReal
numChirps = fileSize/numADCSamples/RX_num/TX_num/frame_num;
LVDS = zeros(1, fileSize);
%create column for each chirp
LVDS = reshape(adcData, numADCSamples*RX_num, numChirps);
%each row is data from one chirp
LVDS = LVDS.';
else
% for complex data
% filesize = 2 * numADCSamples*numChirps
numChirps = fileSize/2/numADCSamples/RX_num/TX_num/frame_num
LVDS = zeros(1, fileSize/2);
%combine real and imaginary part into complex data
%read in file: 2I is followed by 2Q
counter = 1;
for i=1:4:fileSize-1
LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
counter = counter + 2;
end
% create column for each chirp
% LVDS = reshape(LVDS, numADCSamples*RX_num, TX_num*numChirps*frame_num);
% LVDS = reshape(LVDS, numADCSamples, RX_num*numChirps, TX_num, frame_num);
LVDS = reshape(LVDS, numADCSamples, TX_num*RX_num, numChirps, frame_num);
% save LVDS LVDS
%each row is data from one chirp
% LVDS = LVDS.';
% LVDS = permute(LVDS, [2,1,3,4]);
retVal=zeros(RX_num*numChirps,numADCSamples,TX_num, frame_num);
for i = 1:TX_num
    retVal(:,:,i,:)=permute(reshape(LVDS(:,((i-1)*RX_num+1):((i-1)*RX_num+RX_num),:,:), numADCSamples, RX_num*numChirps,frame_num),[2,1,3]);
end
end
%organize data per RX
% adcData = zeros(RX_num,numChirps*numADCSamples*TX_num*frame_num);
% for row = 1: RX_num
%     
%     for f = 1: frame_num
%     
%         for i = 1: numChirps
% 
%             for t = 1: TX_num
% 
%                 adcData(row, (1: numADCSamples) + (t - 1) * numADCSamples + (i - 1) * TX_num * numADCSamples + (f - 1) * TX_num * numChirps * numADCSamples) = LVDS(t + (i - 1) * TX_num + (f - 1) * TX_num * numChirps, (row - 1) * numADCSamples + (1: numADCSamples));
% 
%             end
%             
%         end
%         
%     end
% 
% end

% return receiver data
% retVal = adcData;

end