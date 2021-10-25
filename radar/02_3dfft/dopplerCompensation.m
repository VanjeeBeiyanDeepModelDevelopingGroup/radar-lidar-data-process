%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     dopplerCompensation �����ղ�������                                  %%%
%%%     doppler- doppler FFT�Ľ������Ϊ����������                          %%%
%%%     dopplerIdx- ����������                                              %%%
%%%     TX_num- ������������                                                %%%
%%%     RX_num- ������������                                                %%%
%%%     dopplerBin_num- doppler bin����                                     %%%
%%%     MAX_VEL_ENH_PROCESSING- Max velocityʹ�ܲ���                        %%%
%%%     dopplerCompOut- �����ղ������2D FFT���                            %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.08 version 1.1                            %%%
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