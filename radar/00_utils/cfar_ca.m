%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     CFAR-CA algorithm                                                   %%%
%%%     input- dopplerSum�е�һ�У���һ��range binλ�õ�����doppler bins    %%%                          
%%%     len- doppler bins����                                               %%%
%%%     thresholdScale- ��ֵ                                                %%%
%%%     noiseDivShift- �ο���Ԫ��Ҫ���Ե�ϵ��                               %%%
%%%     guardLen- ������Ԫ����                                              %%%
%%%     winLen- �ο���Ԫ����                                                %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.08 version 1.0                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ index ] = cfar_ca( input, len, thresholdScale, noiseDivShift, guardLen, winLen)
    minInput = min(input);
    maxInput = max(input);
    input = (input-minInput)/(maxInput-minInput);  % normalization
    %initializations
    index = zeros(1, len);
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
        
        total = sumRight + sumLeft;

%         if input(n) > total/power(2, noiseDivShift) + thresholdScale
        if input(n) > total/(2*winLen) + thresholdScale

            index(n) = 1;

        end
        
        n = n + 1;
        
    end
    
    

end