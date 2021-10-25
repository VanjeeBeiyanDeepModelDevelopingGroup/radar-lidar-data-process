function [biMatRAMap,SNRmap] = CFAR2d(RAMap,Tr,Td,Gr,Gd, off_set)
    SNRmap = ones(size(RAMap));
    RDM = RAMap;
    if RDM==zeros(size(RDM))
        biMatRAMap = zeros(size(RDM));
        return
    end
    threshMap = [];
    % *%TODO* :
    %Select the number of Training Cells in both the dimensions.

%     Tr=6;
%     Td=4;

    % *%TODO* :
    %Select the number of Guard Cells in both dimensions around the Cell under 
    %test (CUT) for accurate estimation

%     Gr=2;
%     Gd=2;

    % *%TODO* :
    % offset the threshold by SNR value in dB

%     off_set=0.74;
    % *%TODO* :
    %design a loop such that it slides the CUT across range doppler map by
    %giving margins at the edges for Training and Guard Cells.
    %For every iteration sum the signal level within all the training
    %cells. To sum convert the value from logarithmic to linear using db2pow
    %function. Average the summed values for all of the training
    %cells used. After averaging convert it back to logarithimic using pow2db.
    %Further add the offset to it to determine the threshold. Next, compare the
    %signal under CUT with this threshold. If the CUT level > threshold assign
    %it a value of 1, else equate it to 0.


    % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
    % CFAR

    RDM = RDM/max(max(RDM)); % Normalizing
    [Nr, Nd] = size(RDM);
    % *%TODO* :
    % The process above will generate a thresholded block, which is smaller 
    %than the Range Doppler Map as the CUT cannot be located at the edges of
    %matrix. Hence,few cells will not be thresholded. To keep the map size same
    % set those values to 0. 

    %Slide the cell under test across the complete martix,to note: start point
    %Tr+Td+1 and Td+Gd+1
    for i = Tr+Gr+1:(Nr)-(Tr+Gr)
        for j = Td+Gd+1:(Nd)-(Td+Gd)
            %Create a vector to store noise_level for each iteration on training cells
            noise_level = zeros(1,1);
            %Step through each of bins and the surroundings of the CUT
            for p = i-(Tr+Gr) : i+(Tr+Gr)
                for q = j-(Td+Gd) : j+(Td+Gd)
                    %Exclude the Guard cells and CUT cells
                    if (abs(i-p) > Gr || abs(j-q) > Gd)
                        %Convert db to power
                        noise_level = noise_level + db2pow(RDM(p,q));
                    end
                end
            end

            %Calculate threshould from noise average then add the offset
            threshold = 0;
            if noise_level > 0 || noise_level
                threshold = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
            end
            %Add the SNR to the threshold
            threshold = threshold + off_set;
            threshMap(i,j) = threshold;
            %Measure the signal in Cell Under Test(CUT) and compare against
            CUT = RDM(i,j);

            if (CUT < threshold)
                RDM(i,j) = 0;
            else
                RDM(i,j) = 1;
                SNRmap(i,j) = CUT/threshold;
            end

        end
    end

    RDM(RDM~=0 & RDM~=1) = 0;
    biMatRAMap = RDM;
%     figure;imagesc(threshMap);
end