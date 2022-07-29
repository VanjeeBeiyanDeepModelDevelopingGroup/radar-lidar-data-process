function [outFreqSignal] = musicAlg(signal,L)
%MUSICALG 自适应music算法
%   输入1维信号，按照thresh分割特征向量和特征值，输出频率分析结果
lamda = 1;
d=lamda/2;
[M,~] = size(signal);
% filename1 = ['rawsig_',datestr(now,'HHMMSSFFF'),'.mat'];
% save(filename1,'signal_1d');
% L = 360;
% N:天线个数
N=M;
%% 
% 计算correlation matrix（由于采样长度有限，只能用接收信号的covariance matrix来替代）
% 按照特征向量变化的一阶导来划分信号空间和噪声空间
X = signal*signal';% '这是共轭转置
% filename2 = ['covmat_',datestr(now,'HHMMSSFFF'),'.mat'];
% save(filename2,'X');
[eigVec,eigVal] = eig(X);
Qn = eigVec(:,1:M-1);
%% 扫描所有角度，计算谱分析结果
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
%  所谓steering vector就是目标在某个角度下的对每根天线形成的相位差序列
function [steerV] = array_response_vector(array, theta)
    [N,~] = size(array);
%     # Calculate the steering vector, v for certain angle, theta. Shape of v is N by 1. See equation (5) in [1].
    v = exp(1j*2*pi*array*sin(theta));
%     """ print('The steering vector for certain angle %f is: ' % (theta*180/np.pi), v) """
    steerV = v/sqrt(N);
end