%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     cossp º¯Êý                                                          %%%
%%%                                                                         %%%
%%%     Created by Àî¼Î±¦ 2021.02.08 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Cos = cossp(a)
    
    Zero   =  0.0; 
    MAX    =  1048576.0;
    MIN    =  2.4414062e-4;
    Sign   =  1;
    InvPI  =  0.318309886183791;
    HalfPI =  1.5707963268;
    s4     =  2.601903036e-6;
    s3     = -1.980741872e-4;
    s2     =  8.333025139e-3;
    s1     = -1.666665668e-1;
    C1     =  3.140625;
    C2     =  9.67653589793e-4;

    Y = abs(a) + HalfPI;

    if Y > MAX 
        Y = HalfPI;
    end

    X = Y * InvPI;             %X = Y * (1/PI)         
    N = fix(X);                %N = integer part of X  
    Z = N;                     %Z = float (N)        
  
    if (N&1) ~= 0
        Sign = -Sign;          %quad. 3 or 4
    end

    F = (Y - (Z * C1)) - (Z * C2);      
    R = F;
  
    if F < Zero
        R = -R;
    end

    if R < MIN
        Cos = R * Sign; 
    else
       G = F * F;
       R = (((s4 * G + s3) * G + s2) * G + s1) * G;      
  
       Cos = ((F + F * R) * Sign); 
    end

end