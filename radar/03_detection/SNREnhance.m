function [Xcube_new] = SNREnhance(Xcube,variance)
%SNREnhance �����뵽CFAR֮ǰ�����ݽ��жԱȶ���ǿ
% Xcube_new = 0;
% Xcube = Xcube';
% Xcube = (Xcube-min(Xcube))./(max(Xcube)-min(Xcube));
% Xcube = Xcube/max(Xcube);
%   ���ȷ�����ֱ��ͼ�ֲ���Ȼ��������ӳ����ǿ
% figure(4);
% subplot(4,2,1);imagesc(Xcube);
% subplot(4,2,2);histogram(Xcube);
% %% ����logӳ��
Xcube1 = 20*log10(Xcube);
% subplot(4,2,3);imagesc(Xcube1);
% subplot(4,2,4);histogram(Xcube1);
% title('log ӳ��');
% %% ����gammaӳ��
% gamma = 0.1;
% scale = 1;
% Xcube2 = scale * (Xcube .^ gamma);
% subplot(4,2,5);imagesc(Xcube2);
% subplot(4,2,6);histogram(Xcube2);
% title('gamma ӳ��');
%% ����sigmoidӳ��
miu = max(Xcube1(:))*variance;
alpha = 1/exp(-miu);
k=1:max(Xcube1(:));
y=1./(1+alpha*exp(-k));
Xcube3 = 1./(1+alpha*exp(-Xcube1));
% subplot(4,2,7);imagesc(Xcube3);
% subplot(4,2,8);histogram(Xcube3);
% % subplot(4,2,8);plot(y);
% title('sigmoid ӳ��');
Xcube_new = Xcube3;
%% CLAHE����




end

