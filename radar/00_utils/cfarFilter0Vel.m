function [filteredRows,filteredCols] = cfarFilter0Vel(rows,cols, dopplerBin_num)
%cfarFilter0Vel 仅仅去掉cfar结果里的0速度的index
%   
doppler0idx = dopplerBin_num/2;
filteredRows = [];
filteredCols = [];
for i=1:length(rows)
    if (cols(i)==doppler0idx+2)||(cols(i)==doppler0idx+1)||(cols(i)==doppler0idx)
        continue;
    end
    filteredRows = [filteredRows,rows(i)];
    filteredCols = [filteredCols,cols(i)];
end
end

