function [Xcube_new] = SNREnhance(Xcube,variance)
%SNREnhance 对输入到CFAR之前的数据进行对比度增强
% Xcube_new = 0;
% Xcube = Xcube';
% Xcube = (Xcube-min(Xcube))./(max(Xcube)-min(Xcube));
% Xcube = Xcube/max(Xcube);
%   首先分析其直方图分布，然后对其进行映射增强
% figure(4);
% subplot(4,2,1);imagesc(Xcube);
% subplot(4,2,2);histogram(Xcube);
% %% 进行log映射
Xcube1 = 20*log10(Xcube);
% subplot(4,2,3);imagesc(Xcube1);
% subplot(4,2,4);histogram(Xcube1);
% title('log 映射');
% %% 进行gamma映射
% gamma = 0.1;
% scale = 1;
% Xcube2 = scale * (Xcube .^ gamma);
% subplot(4,2,5);imagesc(Xcube2);
% subplot(4,2,6);histogram(Xcube2);
% title('gamma 映射');
%% 进行sigmoid映射
miu = max(Xcube1(:))*variance;
alpha = 1/exp(-miu);
k=1:max(Xcube1(:));
y=1./(1+alpha*exp(-k));
Xcube3 = 1./(1+alpha*exp(-Xcube1));
% subplot(4,2,7);imagesc(Xcube3);
% subplot(4,2,8);histogram(Xcube3);
% % subplot(4,2,8);plot(y);
% title('sigmoid 映射');
Xcube_new = Xcube3;
%% CLAHE方法




end

