function [fusionConfMap] = fogdetector(path,timeWin)
% function [fusionConfMap] = fogdetector(lidarRAmap,radarRAmap)
% �������۵Ľ������������Ŀ�����˵���ز�ǰ�������Զ�룬��ô�������Ȼ���𽥽��͵�
% ���ǰ���𽥿�������ô����չ����������
% ���Ƕ������������ϵ���˵��������һ����ߵ��ӳ��д��ڶ���޷��ֱ��Ŀ�꣬�����չ��
% ��������д���Ŀ�꣬�ز����ر�Ŀ���������Ϊһ��Ŀ�겻��͸�⣩
% ���������˶���Ŀ�꣬���ں��ײ���˵�����������׼ȷ�����ǽǶȲ���ȷ������Ϊ�Ƿֱ������ޣ�
% ���Ƕ��ڼ����״���˵��һ��Ŀ��ĺ��ر��жϳ�������Ȼû�취��࣬��Ϊǰ����û��չ������ǲ����׼ȷ��
% �����������Ϣ�����ܹ��Ӽ����״�չ���л�ԭ��Ŀ������
fusionConfMap = [];
%% ά��ͳһ����
last_angles_rows = 0;
last_frontNum_cols = 0;
now_angles_rows = 0;
now_frontNum_cols = 0;
%% ��ÿ�������е�5����������ʱ��Ҫ�ó�����
oldest_frontR={};
oldest_delta_frontR={};
oldest_var_frontR={};
oldest_pulseWidth={};
oldest_delta_pulseWidth={};
oldest_var_pulseWidth={};
%% ����ǰ�صı���
last_frontR = [];
new_frontR = [];
sum_frontR = [];
sum_delta_frontR = [];
sum_var_frontR = [];
ave_frontR = [];
ave_delta_frontR = [];
var_frontR = [];
%% ����ı���
last_pulseWidth = [];
new_pulseWidth = [];
sum_pulseWidth = [];
sum_delta_pulseWidth = [];
sum_var_pulseWidth = [];
ave_pulseWidth = [];
ave_delta_pulseWidth = [];
var_pulseWidth = [];
%% ������
timeWinCnt = 0;
outputCnt = 0;
for i=1:100 % �ӵڶ�֡��ʼ
    %% load����
    [lidarRAmap,radarRAmap] = loadRAmaps(i,path);
    [rangeIdxFrontR,clusterPulseWidth] = lidarFnWcal(lidarRAmap,15,10,0.25);
    %% ������֡������ά��ͳһ�������ۼ�����ά��ҲҪͳһ������д��̫�鷳����matlab�����鷳
    [now_angles_rows,now_frontNum_cols] = size(rangeIdxFrontR);
    % ��ǰ֡����һ֡ά�ȴ�
    if now_frontNum_cols>last_frontNum_cols
        % ����ά��֮��Ҫ����ά�ȴ�С��������ǲ���Ҫ�ģ���Ϊ����С�Ĳ������ȥ��
%         now_frontNum_cols = last_frontNum_cols;
        % ǰ��
        last_frontR_fill = 564*ones(now_angles_rows,now_frontNum_cols);
        last_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = last_frontR;
        last_frontR = last_frontR_fill;
        % ����
        last_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        last_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = last_pulseWidth;
        last_pulseWidth = last_pulseWidth_fill;
        % ǰ��λ���ۼ�����sum_frontR
        sum_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_frontR;
        sum_frontR = sum_frontR_fill;
        % ǰ�ر仯���ۼ�����sum_delta_frontR
        sum_delta_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_delta_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_delta_frontR;
        sum_delta_frontR = sum_delta_frontR_fill;
        % ǰ�ط����ۼ�����sum_var_frontR
        sum_var_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_var_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_var_frontR;
        sum_var_frontR = sum_var_frontR_fill;
        % �����С�ۼ�����sum_pulseWidth
        sum_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_pulseWidth;
        sum_pulseWidth = sum_pulseWidth_fill;
        % ����仯���ۼ�����sum_delta_pulseWidth
        sum_delta_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_delta_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_delta_pulseWidth;
        sum_delta_pulseWidth = sum_delta_pulseWidth_fill;
        % �������ۼ�����sum_var_pulseWidth
        sum_var_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        sum_var_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = sum_var_pulseWidth;
        sum_var_pulseWidth = sum_var_pulseWidth_fill;
        % ǰ��λ�þ�ֵ
        ave_frontR_fill = zeros(now_angles_rows,now_frontNum_cols);
        ave_frontR_fill(1:last_angles_rows,1:last_frontNum_cols) = ave_frontR;
        ave_frontR = ave_frontR_fill;
        % �����С��ֵ
        ave_pulseWidth_fill = zeros(now_angles_rows,now_frontNum_cols);
        ave_pulseWidth_fill(1:last_angles_rows,1:last_frontNum_cols) = ave_pulseWidth;
        ave_pulseWidth = ave_pulseWidth_fill;
    % ��ǰ֡����һ֡ά��С
    elseif now_frontNum_cols<last_frontNum_cols
        % ǰ��
        rangeIdxFrontR_fill = 564*ones(last_angles_rows,last_frontNum_cols);
        rangeIdxFrontR_fill(1:now_angles_rows,1:now_frontNum_cols) = rangeIdxFrontR;
        rangeIdxFrontR = rangeIdxFrontR_fill;
        % ����
        clusterPulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
        clusterPulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = clusterPulseWidth;
        clusterPulseWidth = clusterPulseWidth_fill;
%         % ǰ��λ���ۼ�����sum_frontR
%         sum_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_frontR;
%         sum_frontR = sum_frontR_fill;
%         % ǰ�ر仯���ۼ�����sum_delta_frontR
%         sum_delta_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_delta_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_delta_frontR;
%         sum_delta_frontR = sum_delta_frontR_fill;
%         % ǰ�ط����ۼ�����sum_var_frontR
%         sum_var_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_var_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_var_frontR;
%         sum_var_frontR = sum_var_frontR_fill;
%         % �����С�ۼ�����sum_pulseWidth
%         sum_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_pulseWidth;
%         sum_pulseWidth = sum_pulseWidth_fill;
%         % ����仯���ۼ�����sum_delta_pulseWidth
%         sum_delta_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_delta_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_delta_pulseWidth;
%         sum_delta_pulseWidth = sum_delta_pulseWidth_fill;
%         % �������ۼ�����sum_var_pulseWidth
%         sum_var_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         sum_var_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = sum_var_pulseWidth;
%         sum_var_pulseWidth = sum_var_pulseWidth_fill;
%         % ǰ��λ�þ�ֵ
%         ave_frontR_fill = zeros(last_angles_rows,last_frontNum_cols);
%         ave_frontR_fill(1:now_angles_rows,1:now_frontNum_cols) = ave_frontR;
%         ave_frontR = ave_frontR_fill;
%         % �����С��ֵ
%         ave_pulseWidth_fill = zeros(last_angles_rows,last_frontNum_cols);
%         ave_pulseWidth_fill(1:now_angles_rows,1:now_frontNum_cols) = ave_pulseWidth;
%         ave_pulseWidth = ave_pulseWidth_fill;
        % ����ά��֮��Ҫ����ά�ȴ�С��������ǲ���Ҫ�ģ���Ϊ����С�Ĳ������ȥ��
        now_frontNum_cols = last_frontNum_cols;
    % �����ͬ��ֱ��������
    end
    %% ǰ��
    new_frontR = rangeIdxFrontR;
    sum_frontR = sum_frontR+new_frontR;
    sum_delta_frontR = sum_delta_frontR+(new_frontR-last_frontR);
    sum_var_frontR = sum_var_frontR+abs(new_frontR-ave_frontR);
    % ����
    oldest_frontR{end+1} = new_frontR;
    oldest_delta_frontR{end+1} = new_frontR-last_frontR;
    oldest_var_frontR{end+1} = abs(new_frontR-ave_frontR);
    %% ����
    new_pulseWidth = clusterPulseWidth;
    sum_pulseWidth = sum_pulseWidth+new_pulseWidth;
    sum_delta_pulseWidth = sum_delta_pulseWidth+(new_pulseWidth-last_pulseWidth);
    sum_var_pulseWidth = sum_var_pulseWidth+abs(new_pulseWidth-ave_pulseWidth);
    % ����
    oldest_pulseWidth{end+1} = new_pulseWidth;
    oldest_delta_pulseWidth{end+1} = new_pulseWidth-last_pulseWidth;
    oldest_var_pulseWidth{end+1} = abs(new_pulseWidth-ave_pulseWidth);
    %% ��������û���壩
%     timeWinCnt = timeWinCnt+1;
    if i > timeWin
        %% ����ؼ���������
        ave_frontR = sum_frontR/timeWin; % ǰ�ؾ�ֵλ��
        ave_delta_frontR = sum_delta_frontR/timeWin; % ǰ�ش��ھ�ֵ�仯��
        var_frontR = sum_var_frontR/timeWin; % ǰ�ط���
        
        ave_pulseWidth = sum_pulseWidth/timeWin; % �����ֵ
        ave_delta_pulseWidth = sum_delta_pulseWidth/timeWin;% ����仯��
        var_pulseWidth = sum_var_pulseWidth/timeWin; % ������
        
        %% �����ۼӱ���
        sum_frontR  = sum_frontR-oldest_frontR{1};
        sum_delta_frontR = sum_delta_frontR-oldest_delta_frontR{1};
        sum_var_frontR = sum_var_frontR-oldest_var_frontR{1};
        
        sum_pulseWidth = sum_pulseWidth-oldest_pulseWidth{1};
        sum_delta_pulseWidth = sum_delta_pulseWidth-oldest_delta_pulseWidth{1};
        sum_var_pulseWidth = sum_var_pulseWidth-oldest_var_pulseWidth{1};
        
        %% ��յ����ϵ��Ǹ�
        oldest_frontR(1)=[];
        oldest_delta_frontR(1)=[];
        oldest_var_frontR(1)=[];
        oldest_pulseWidth(1)=[];
        oldest_delta_pulseWidth(1)=[];
        oldest_var_pulseWidth(1)=[];
    end
    %% ����
    last_frontR = new_frontR;
    last_pulseWidth = new_pulseWidth;
    last_angles_rows=now_angles_rows;
    last_frontNum_cols=now_frontNum_cols;
    %% ���������ж�
    % ave_frontR = sum_frontR/timeWinCnt; % ǰ�ؾ�ֵλ��
    % ave_delta_frontR = sum_delta_frontR/timeWinCnt; % ǰ�ش��ھ�ֵ�仯��
    % var_frontR = sum_var_frontR/timeWinCnt; % ǰ�ط���        
    % ave_pulseWidth = sum_pulseWidth/timeWinCnt; % �����ֵ
    % ave_delta_pulseWidth = sum_delta_pulseWidth/timeWinCnt;% ����仯��
    % var_pulseWidth = sum_var_pulseWidth/timeWinCnt; % ������
    if isempty(ave_frontR) || isempty(ave_delta_frontR) || isempty(var_frontR) || isempty(ave_pulseWidth) || isempty(ave_delta_pulseWidth) || isempty(var_pulseWidth)
        fprintf("���ŵ�%d֡�����ڻ�û�����\n",i);
        continue
    end
    [angleRow,peakNum] = size(ave_frontR); % ���ĸ���һ����ά�ȶ���һ����
%     for angle=1:angleRow
%         for num=1:peakNum
%             %% ǰ��
    
    
end
end

