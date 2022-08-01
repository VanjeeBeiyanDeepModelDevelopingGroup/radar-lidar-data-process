function [outFreqSignal] = musicAlg(signal,L)
%MUSICALG ����Ӧmusic�㷨
%   ����1ά�źţ�����thresh�ָ���������������ֵ�����Ƶ�ʷ������
lamda = 1;
d=lamda/2;
[M,~] = size(signal);
% filename1 = ['rawsig_',datestr(now,'HHMMSSFFF'),'.mat'];
% save(filename1,'signal_1d');
% L = 360;
% N:���߸���
N=M;
%% 
% ����correlation matrix�����ڲ����������ޣ�ֻ���ý����źŵ�covariance matrix�������
% �������������仯��һ�׵��������źſռ�������ռ�
X = signal*signal';% '���ǹ���ת��
% filename2 = ['covmat_',datestr(now,'HHMMSSFFF'),'.mat'];
% save(filename2,'X');
[eigVec,eigVal] = eig(X);
Qn = eigVec(:,1:M-1);
%% ɨ�����нǶȣ������׷������
antennaArr = linspace(0,(N-1)*d,N)';
% scanAngle = linspace(-pi/2,pi/2,L)';
scanAngle = linspace(pi/2,-pi/2,L)';
powerSpectrumInSpace = zeros(1,L);
for i =1:L
    av = array_response_vector(antennaArr,scanAngle(i));
    powerSpectrumInSpace(i) = 1/norm(Qn'*av);
end
outFreqSignal = (powerSpectrumInSpace-min(powerSpectrumInSpace))/(max(powerSpectrumInSpace)-min(powerSpectrumInSpace));
% outFreqSignal = powerSpectrumInSpace/max(powerSpectrumInSpace);

end
%  Calculate the steering vector for a certain angle, theta.
%  ��νsteering vector����Ŀ����ĳ���Ƕ��µĶ�ÿ�������γɵ���λ������
function [steerV] = array_response_vector(array, theta)
    [N,~] = size(array);
%     # Calculate the steering vector, v for certain angle, theta. Shape of v is N by 1. See equation (5) in [1].
    v = exp(1j*2*pi*array*sin(theta));
%     """ print('The steering vector for certain angle %f is: ' % (theta*180/np.pi), v) """
    steerV = v/sqrt(N);
end