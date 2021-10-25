%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     dopplerCompensation 多普勒补偿函数                                  %%%
%%%     doppler- doppler FFT的结果，作为函数的输入                          %%%
%%%     dopplerIdx- 多普勒索引                                              %%%
%%%     TX_num- 发射天线数量                                                %%%
%%%     RX_num- 接收天线数量                                                %%%
%%%     dopplerBin_num- doppler bin数量                                     %%%
%%%     MAX_VEL_ENH_PROCESSING- Max velocity使能参数                        %%%
%%%     dopplerCompOut- 多普勒补偿后的2D FFT结果                            %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.02.08 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function doppler = dopplerCompensation( doppler, dopplerIdx, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING)
    
    dopplerSignIdx = dopplerIdx - (dopplerBin_num / 2);

    if(MAX_VEL_ENH_PROCESSING)
        numHypotheses = TX_num;
        wrapInd = - fix(numHypotheses / 2);
            
        if (TX_num == 2) && (dopplerSignIdx < 0)
        % e.g.2 tx, doppler < 0: wrapInd = 0 
        
            wrapInd = wrapInd + 1;
            
        end
    else
        
            numHypotheses = 1;
            wrapInd = 0;
    end

    numVirtualAnt = TX_num * RX_num;

    for hypothesisIdx = 0: (numHypotheses - 1)
  
        virtAntIdx = 1;
        % transfer data corresponding to azimuth virtual antennas (corresponding to chirp of antenna Tx0)
        for j = 1: RX_num
        
            doppler(virtAntIdx + hypothesisIdx * numVirtualAnt) = doppler(virtAntIdx);
            virtAntIdx = virtAntIdx + 1;
        end

        if TX_num > 1
        
            dopplerCompensationIdx = (dopplerSignIdx + (wrapInd * dopplerBin_num)) / TX_num;
            Cos = cossp(2 * pi * dopplerCompensationIdx / dopplerBin_num);
            Sin = sinsp(2 * pi * dopplerCompensationIdx / dopplerBin_num);

            for txAntIdx = 1: TX_num - 1
            
                for j = 1: RX_num
                
                    doppler(virtAntIdx + hypothesisIdx * numVirtualAnt) = AoAProcDSP_dopplerComp(doppler(virtAntIdx), Cos, Sin);
                    virtAntIdx = virtAntIdx + 1;
                    
                end
                if txAntIdx < (TX_num - 1)
                
                    %Increment Doppler phase shift
                    temp = Cos * Cos - Sin * Sin;
                    Sin = 2 * Cos * Sin;
                    Cos = temp;
                    
                end
                
            end
        
        end

        wrapInd = wrapInd + 1;
        
    end

end