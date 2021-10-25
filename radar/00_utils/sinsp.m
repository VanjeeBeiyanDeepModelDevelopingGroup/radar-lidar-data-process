%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     sinsp ����                                                          %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.08 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Sin = sinsp(a)
    
    InvPI =  0.318309886183791;
    One   =  1.0;
    MAX   =  1048576.0;
    Zero  =  0.0;
    s1    = -1.666665668e-1;
    s2    =  8.333025139e-3; 
    s3    = -1.980741872e-4;
    s4    =  2.601903036e-6; 
    C1    =  3.140625;
    C2    =  9.67653589793e-4;

    Sign = One;
    Y    = a;

    if abs(Y) > MAX 
        Y = Zero;
    end

    X = Y * InvPI;         % X = Y * (1/PI) 
    N = fix(X);            % N = integer part of X 
    Z = N;                         
  
    if (N & 1) ~= 0
        Sign = - Sign;           % Quadrant 3 or 4 
    end

    F = (Y - Z * C1) - Z * C2;      
    G = F * F;
    R = (((s4 * G + s3) * G + s2) * G + s1) * G;
        
    Sin = (F + F * R) * Sign;

end