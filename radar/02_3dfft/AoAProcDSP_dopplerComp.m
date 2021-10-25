%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AoAProcDSP_dopplerComp ����                                         %%%
%%%     input- ����������                                                   %%%
%%%     Cos, Sin- ������ϵ��                                                %%%
%%%     Output- �����Ľ�����                                              %%%
%%%                                                                         %%%
%%%     Created by ��α� 2021.02.08 version 1.1                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Output = AoAProcDSP_dopplerComp(input, Cos, Sin)

    yRe = real(input) * Cos + imag(input) * Sin;
    yIm = imag(input) * Cos - real(input) * Sin;
    Output = yRe + 1i * yIm;

end