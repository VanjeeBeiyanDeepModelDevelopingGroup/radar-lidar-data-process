%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AoAProcDSP_dopplerComp 函数                                         %%%
%%%     input- 函数的输入                                                   %%%
%%%     Cos, Sin- 正余弦系数                                                %%%
%%%     Output- 函数的结果输出                                              %%%
%%%                                                                         %%%
%%%     Created by 李嘉宝 2021.02.08 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = AoAProcDSP_dopplerComp(input, Cos, Sin)

    yRe = real(input) * Cos + imag(input) * Sin;
    yIm = imag(input) * Cos - real(input) * Sin;
    Output = yRe + 1i * yIm;

end