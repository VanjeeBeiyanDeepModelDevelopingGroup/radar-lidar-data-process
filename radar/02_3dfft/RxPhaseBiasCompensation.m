%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     RxPhaseBiasCompensation 接收天线相位补偿函数                                                     %%%
%%%     angleIn- 根据CFAR结果选取的虚拟天线个数doppler FFT的结果值，作为函数的输入                       %%%                          
%%%     TX_num- 发射天线数量                                                                             %%%
%%%     RX_num- 接收天线数量                                                                             %%%
%%%     MAX_VEL_ENH_PROCESSING- max velocity是否使能标志                                                 %%%
%%%     AOAIn- 接收天线相位补偿后的值作为AOA FFT的输入                                                   %%%
%%%                                                                                                      %%%
%%%     Created by 李嘉宝 2021.02.08 version 1.0                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AOAIn = RxPhaseBiasCompensation( angleIn, TX_num, RX_num, MAX_VEL_ENH_PROCESSING)
    AOAIn = zeros(1, TX_num * RX_num);

    if MAX_VEL_ENH_PROCESSING
        
        rxCoef = [5863, 19158; 2858, 23451; 4568, 22093; 0, 16399];
        
    else
        
        rxCoef = [21646, 7035; 21172, 18420; 15490, 21080; -3905, 25130; 0, 16399; -7985, 18443; -10962, 15291; -17653, 3133; 386, -17208; 8587, -18744; 11857, -15772; 18493, -2907];
        
    end
       
    for m = 1: TX_num

        for n = 1: RX_num

            Real = round(rxCoef(n + (m - 1) * RX_num, 2) * real(angleIn(n + (m - 1) * RX_num))) - round(rxCoef(n + (m - 1) * RX_num, 1) * imag(angleIn(n + (m - 1) * RX_num)));
            Imag = round(rxCoef(n + (m - 1) * RX_num, 2) * imag(angleIn(n + (m - 1) * RX_num))) + round(rxCoef(n + (m - 1) * RX_num, 1) * real(angleIn(n + (m - 1) * RX_num)));
            AOAIn(n + (m - 1) * RX_num) = Real + Imag * 1i;

        end

    end

end