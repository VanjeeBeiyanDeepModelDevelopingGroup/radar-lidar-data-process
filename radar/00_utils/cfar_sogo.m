%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     CFAR-SOGO algorithm                                                 %%%
%%%     input- dopplerSum中的一行，即一个range bin位置的所有doppler bins    %%%                          
%%%     len- range bins数量                                                 %%%
%%%     thresholdScale- 阈值                                                %%%
%%%     noiseDivShift- 参考单元需要除以的系数                               %%%
%%%     guardLen- 保护单元数量                                              %%%
%%%     winLen- 参考单元数量                                                %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.02.08 version 1.0                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ index ] = cfar_sogo( input, len, cfartype, thresholdScale, noiseDivShift, guardLen, winLen)

    %initializations
    index = zeros(len, 1);
    n = 1;
    
    sumLeft = sum(input(len - guardLen - (0: winLen - 1)));
    sumRight = sum(input(guardLen + (2: winLen + 1)));
    
    leftIdx = len - guardLen;
    rightIdx = 1 + guardLen + winLen;
    
    while n <= len
        
        if n > 1
           
            if leftIdx <= len
                
                sumLeft = sumLeft + input(leftIdx) - input(leftIdx - winLen);
            
            elseif leftIdx <= (len + winLen)
                
                sumLeft = sumLeft + input(leftIdx - len) - input(leftIdx - winLen);
                
                if leftIdx == len + winLen
                    
                    leftIdx = winLen;
                    
                end
                
            end
            
        end
        
        leftIdx = leftIdx + 1;
        
        if n > 1
            
            if rightIdx <= len
                
                sumRight = sumRight + input(rightIdx) - input(rightIdx - winLen);
                
            elseif rightIdx <= (len + winLen)
                
                sumRight = sumRight + input(rightIdx - len) - input(rightIdx - winLen);
                
                if rightIdx == len + winLen
                    
                    rightIdx = winLen;
                    
                end
                
            end
            
        end
        
        rightIdx = rightIdx + 1;
        
        if n <= guardLen + winLen
            
            if input(n) > sumRight/power(2, noiseDivShift - 1) + thresholdScale
                
                index(n) = 1;
                
            end
            
        elseif n > len - (guardLen + winLen)
            
            if input(n) > sumLeft/power(2, noiseDivShift - 1) + thresholdScale
                
                index(n) = 1;
                
            end
            
        else
            
            if strcmp(cfartype, 'ca')
                
                if input(n) > (sumLeft + sumRight)/power(2, noiseDivShift) + thresholdScale
                
                    index(n) = 1;
                
                end
                
            elseif strcmp(cfartype, 'so')
                
                if sumLeft < sumRight
                    
                    if input(n) > sumLeft/power(2, noiseDivShift - 1) + thresholdScale
                
                        index(n) = 1;
                
                    end
                    
                else
                    
                    if input(n) > sumRight/power(2, noiseDivShift - 1) + thresholdScale
                
                        index(n) = 1;
                
                    end
                    
                end
                
            elseif strcmp(cfartype, 'go')
                
                if sumLeft > sumRight
                    
                    if input(n) > sumLeft/power(2, noiseDivShift - 1) + thresholdScale
                
                        index(n) = 1;
                
                    end
                    
                else
                    
                    if input(n) > sumRight/power(2, noiseDivShift - 1) + thresholdScale
                
                        index(n) = 1;
                
                    end
                    
                end
                
            end
            
        end
        
        n = n + 1;
        
    end

end