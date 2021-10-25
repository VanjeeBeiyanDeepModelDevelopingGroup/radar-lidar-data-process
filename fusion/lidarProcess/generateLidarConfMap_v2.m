% ���µ����������µļ��Ⲩ�����ŷֲ�ͼ
% ��ʵ�������������ǣ���ȡ->��������->��������
function [lidarConfMap,lidarMiuSig] = generateLidarConfMap_v2(LidarADmap,lidarAngleGrid,lidarRangeGrid)

[angles,ranges] = size(LidarADmap);
lidarConfMap = zeros(angles,ranges);
lidarMiuSig = [];
for j=1:angles
    singleLidarAD = LidarADmap(j,:);
    %% ȡ����
    bit_singleLidarAD = cfar_ca1D_square(singleLidarAD,50,25,0.05,0);
%     subplot(2,3,[2,3]);plot(50*bit_singleLidarAD+130);
    len = length(bit_singleLidarAD);
    old_bit_ = bit_singleLidarAD(1);
    old_val = singleLidarAD(1);
    err_val = zeros(1,len);
    peakStartEnd = [];
    startF = 1;
    endF = 1;
    peakValue = 0;
    peakIdx = 1;
    peakArray = [];
    for range=2:len
        % ��cfar��һ�׵�
        new_bit_ = bit_singleLidarAD(range);
        delta_bit_ = new_bit_-old_bit_;
        old_bit_ = new_bit_;
%         % ��ԭʼ�źŵ�һ�׵�
%         new_val = singleLidarAD(range);
%         delta_val(range) = new_val-old_val;
%         old_val = new_val;
        % ������λ��
        if delta_bit_ == 1
            startF = range;
            % ��ֵ��С��λ��
            value = singleLidarAD(range);
            if value > peakValue
                peakValue = value;
                peakIdx = range;
            end
        % �½���λ��
        elseif delta_bit_ == -1
            endF = range;
            peakStartEnd = [peakStartEnd,[startF,endF]];
            startF = 1;
            endF = 1;
            % �����½���֮��Ѽ�ֵ��С������������
            peakArray = [peakArray,[peakValue,peakIdx]];
            peakValue = 0;
            peakIdx = 1;
        end
        value = singleLidarAD(range);
        if value > peakValue
            peakValue = value;
            peakIdx = range;
        end
    end
    crestNum = length(peakStartEnd)/2;
    % �ڲ����б���
    for crest=1:crestNum
        startF = peakStartEnd(2*crest-1);
        endF = peakStartEnd(2*crest);
        crestArray = singleLidarAD(startF:endF);
        crestCnts = length(crestArray);
        old_val = crestArray(1);
        old_error_val = 0;
        err_err_val = zeros(1,length(singleLidarAD));
        err_val = zeros(1,length(singleLidarAD));
        peakVal = 0;
        peakIdx = 1;
        for crestCnt = 2:crestCnts
            %% �󲨷��һ�׵�����׵�
            new_val = crestArray(crestCnt);
            new_error_val = new_val-old_val;
            err_val(startF+crestCnt-2) = new_error_val;
            err_err_val(startF+crestCnt-2) = new_error_val-old_error_val;
            %% �󲨷��һ�׵���ֵ����λ�ã�Ҳ�������еĶ��׵����λ�ã����׵������ɢ̫�����ˣ�ȡ����������Ȼ��Сһ�׵���ֵ�Ĳ�Ҫ
            if new_error_val > peakVal
                peakVal = new_error_val;
                peakIdx = startF+crestCnt-2;
            end
            %% ����
            old_val = new_val;
            old_error_val = new_error_val;
        end
        %% ����һ�׵���ֵ������ſ�ʼ��λ�õı���Ϊǰ�ط���
%         frontScore = (peakIdx+2-startF)/crestCnts;
        startRange = lidarRangeGrid(startF)/lidarRangeGrid(end);
%         score = peakVal//(endF-peakIdx-2); % ��ֱ�ӵ�����
%         score = peakVal*startRange^2/(endF-peakIdx-2); % ��ò���ѽ
        score = peakVal*startRange/(peakIdx+2-startF);
        if score == 0
            continue;
        end
        sigma = 0.5/score;
%         sigma = 30*atan(score);
%         if endF-peakIdx-2 >0
        if peakIdx+2-startF > 0
            angle0 = lidarAngleGrid(j);
            range0 = lidarRangeGrid(peakIdx);
            lidarMiuSig = [angle0,range0,sigma;lidarMiuSig];% ����angle,range,score
        %% ���㼤��Ŀ���ֲ�ͼ
            for m=1:angles
                for n=1:ranges
%                     sigma = sigma;
                    distance = (((j-m))^2+(peakIdx-n)^2)/(sigma^2);
                    if distance > 36 % �޶��ֲ���Χ
                        distance = 36;
                    end
                    value = exp(-distance/2);
                    % ��������ֵ����ԭ����ֵ���滻��ԭ����
                    if value > lidarConfMap(m,n)
                        lidarConfMap(m,n) = value;
                    end
                end
            end
        end
    end
end

end

