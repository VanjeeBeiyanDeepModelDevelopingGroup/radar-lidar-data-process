%%% Created by 李嘉宝 2021.05.06 version 1.2 %%%

function read_data2_v1_4
clear
close all
f = figure('Visible','off','Position',[30,50,1180,600]);
plot_data = axes('Units','Pixels','Position',[50,300,800,200]);
plot_mmwave_data = axes('Units','Pixels','Position',[50,50,800,200]);
%  Construct the components.
hfile = uicontrol('Style','pushbutton','Units','pixels',...
    'String','激光文件','FontSize',12,...
    'Position',[50,555,100,30],...
    'Callback',@file_button_Callback);
% htext = uicontrol(f,'Style','text','Units','pixels',...
%     'FontSize',12,'Position',[180,560,200,30],...
%     'Background','w','HorizontalAlignment','left');
% rangetext = uicontrol(f,'Style','text','Units','pixels',...
%     'FontSize',12,'Position',[430,560,200,30],...
%     'Background','w','HorizontalAlignment','left');
frameNextPackButton = uicontrol('Style','pushbutton','Units','pixels',...
    'String','next','FontSize',12,...
    'Position',[180,555,40,30],...
    'Callback',@frame_change_next_Callback);
frameResetButton = uicontrol('Style','pushbutton','Units','pixels',...
    'String','reset','FontSize',12,...
    'Position',[230,555,40,30],...
    'Callback',@frame_change_reset_Callback);    
mmwavefile = uicontrol('Style','pushbutton','Units','pixels',...
    'String','毫米波文件','FontSize',12,...
    'Position',[50,516,100,30],...
    'Callback',@mmwavefile_Callback);
frametext = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[177,510,130,30],...
    'String','输入开始索引：','HorizontalAlignment','left');
frameIndex = uicontrol('Style','edit','Units','pixels',...
    'String','','FontSize',12,...
    'Position',[300,515,35,30]);
% frameNumtext = uicontrol(f,'Style','text','Units','pixels',...
%     'FontSize',12,'Position',[668,555,200,30],...
%     'Visible','off','HorizontalAlignment','left');
videoFile = uicontrol('Style','pushbutton','Units','pixels',...
    'String','视频文件','FontSize',12,...
    'Position',[50,500,100,30],...
    'Callback',@video_button_Callback);
framePlotButton = uicontrol('Style','pushbutton','Units','pixels',...
    'String','start','FontSize',12,...
    'Position',[350,516,40,30],...
    'Callback',@frame_change_disp_Callback);
frameStopButton = uicontrol('Style','pushbutton','Units','pixels',...
    'String','stop','FontSize',12,...
    'Position',[395,516,40,30],...
    'Callback',@frame_change_stop_Callback);  
modeButtonGroup = uibuttongroup('Visible','on','Units','pixels',...
                  'Position',[500 516 180 40],'Title','Modes',...
                  'SelectionChangedFcn',@bselection);
mode1 = uicontrol(modeButtonGroup,'Style','radiobutton',...
                  'String','urr',...
                  'Position',[10 10 30 15],...
                  'HandleVisibility','on');
mode2 = uicontrol(modeButtonGroup,'Style','radiobutton',...
                  'String','mrr',...
                  'Position',[60 10 50 15],...
                  'HandleVisibility','on'); 
mode3 = uicontrol(modeButtonGroup,'Style','radiobutton',...
                  'String','mrr+urr',...
                  'Position',[110 10 70 15],...
                  'HandleVisibility','on');
methodButtonGroup1 = uibuttongroup('Visible','off','Units','pixels',...
                  'Position',[700 516 105 80],'Title','Methods');
method1 = uicontrol(methodButtonGroup1,'Style','radiobutton',...
                  'String','test',...
                  'Position',[10 50 50 15]);
method2 = uicontrol(methodButtonGroup1,'Style','radiobutton',...
                  'String','doppler comp',...
                  'Position',[10 30 100 15]);
method3 = uicontrol(methodButtonGroup1,'Style','radiobutton',...
                  'String','tracking',...
                  'Position',[10 10 70 15]); 
methodButtonGroup2 = uibuttongroup('Visible','off','Units','pixels',...
                  'Position',[700 516 105 80],'Title','Methods');
method1 = uicontrol(methodButtonGroup2,'Style','radiobutton',...
                  'String','normal',...
                  'Position',[10 50 50 15]);
method2 = uicontrol(methodButtonGroup2,'Style','radiobutton',...
                  'String','CRT',...
                  'Position',[10 30 100 15]);
method3 = uicontrol(methodButtonGroup2,'Style','radiobutton',...
                  'String','apply',...
                  'Position',[10 10 70 15]);              
finalResultTable = uitable(f, 'ColumnName', {'Distance'; 'Velocity'; 'x'; 'y'; 'z'}, 'Position', [880, 50, 275, 500]);

if strcmp(modeButtonGroup.SelectedObject.String, 'urr')
    methodButtonGroup1.Visible = 'on';
end

hfile.Units = 'normalized';
plot_data.Units = 'normalized';
plot_mmwave_data.Units = 'normalized';
htext.Units = 'normalized';
rangetext.Units = 'normalized';
mmwavefile.Units = 'normalized';
frametext.Units = 'normalized';
frameIndex.Units = 'normalized';
frameNumtext.Units = 'normalized';
framePlotButton.Units = 'normalized';
frameStopButton.Units = 'normalized';
frameNextPackButton.Units = 'normalized';
frameResetButton.Units = 'normalized';
finalResultTable.Units = 'normalized';
modeButtonGroup.Units = 'normalized';
methodButtonGroup.Units = 'normalized';
f.Units = 'normalized';
f.Name = '数据显示_version:4';
movegui(f,'center')
f.Visible = 'on';
global lidar_dir_path
global radar_dir_path
global video_dir_path
global radarDataAll
global vel1data
global vel2data
global mmwavefile_name
% global lidarData;
global lidarDataFrame
global filename_pack
global nextCnt
global frame_num
global stopFlag
global area_cor_form4
global x_cfd_cor_form4




    function video_button_Callback(~,~)
        %% 读取视频文件路径
        [video_dir_path] = uigetdir()
    end

    function file_button_Callback(~,~)
        [lidar_dir_path] = uigetdir()
%         fprintf("开始读激光雷达数据......\n");
%         filename
%         data = load([dir_path,filename]);
%         lidarData = data.data;
% %         lidarDataFrame = packLidarDataIntoFrames(lidarData); % 角度*距离*帧号
%         lidarDataFrame = packLidarDataIntoFrames_v02(lidarData); % 角度*距离*帧号
%         size(lidarDataFrame)
%         fprintf("结束!\n");
%         fprintf("读取修正参数......\n");
        %% 载入修正参数
%         area_t = load('.\lidar\02_detection\rangeMeasure\measureParam\4_area_t.mat');%载入通道4修正系数
        area_t = load('./lidar/02_detection/rangeMeasure/measureParam/4_area_t.mat');%载入通道4修正系数,unix系统中是反向斜线
        area_cor_form4=area_t.area_t(1,:);
        x_cfd_cor_form4=area_t.area_t(2,:);
        % load('.\measureParam\5_area_t.mat');%载入通道5修正系数
        % area_cor_form5=area_t(1,:);
        % x_cfd_cor_form5=area_t(2,:);
        % load('.\measureParam\6_area_t.mat');%载入通道6修正系数
        % area_cor_form6=area_t(1,:);
        % x_cfd_cor_form6=area_t(2,:);
        %% 直接对其他路径进行赋值
        video_dir_path = lidar_dir_path(1:end-5)
        radar_dir_path = [video_dir_path,'radar']
    end
    %% 读取毫米波文件
    function mmwavefile_Callback(~,~)         
     radar_dir_path = uigetdir()
    end
    
    function bselection(~,event)
        if strcmp(event.NewValue.String, 'urr') || strcmp(modeButtonGroup.SelectedObject.String, 'mrr')
            methodButtonGroup2.Visible = 'off';
            methodButtonGroup1.Visible = 'on';
        elseif strcmp(event.NewValue.String, 'mrr+urr')
            methodButtonGroup1.Visible = 'off';
            methodButtonGroup2.Visible = 'on';
        end
    end

    function frame_change_disp_Callback(~,~)
        lidar_sort_nat_name = dir(lidar_dir_path);
        radar_sort_nat_name = dir(radar_dir_path);
        lidar_dir_struct=sort_nat({lidar_sort_nat_name.name});
        radar_dir_struct=sort_nat({radar_sort_nat_name.name});
        frameNum = length(lidar_dir_struct);
        timestampShift = 5;
        startFrameNum = str2num(frameIndex.String);
        % 获取视频文件
        vid = VideoReader([video_dir_path, '/video.mp4']);
        if vid == -1
           fprintf('不存在对应的视频文件！\n');
           return;
        end
        stopFlag = false;
        for i=startFrameNum+2:frameNum
            if stopFlag
                fprintf("停止!\n");
                break;
            end
            fprintf('帧号：%d\n',i);
        %% 获取文件名
%             lidar_filename = [lidar_dir_path,'\', lidar_dir_struct{i}];
%             radar_filename = [radar_dir_path,'\', radar_dir_struct{i}];
            lidar_filename = [lidar_dir_path,'/', lidar_dir_struct{i}]; % unix系统中是反向斜线
            radar_filename = [radar_dir_path,'/', radar_dir_struct{i}];% unix系统中是反向斜线
        %% 将帧号对齐
%             lidarFrameNum = str2num(lidar_dir_struct{i}(1:end-4));
%             radarFrameNum = str2num(radar_dir_struct{i}(1:end-4));
%             frameDiff = lidarFrameNum-radarFrameNum;
%             checkCnt = 0;
%             if frameDiff~=0
%                 % 向后(或前)查看timestampShift帧
%                 for checkCnt = 1:timestampShift
%                     if frameDiff < 0 && frameDiff > timestampShift
%                         checkCnt = checkCnt*(-1);
%                     end
%                     radarFrameNum_temp = str2num(radar_dir_struct{i+checkCnt}(1:end-4));
%                     % 找到了完全相同的帧号，就用这个文件
%                     if lidarFrameNum == radarFrameNum_temp
%                         radar_filename = [radar_dir_path,'\', radar_dir_struct{i+checkCnt}];
%                         break
%                     end
%                 end
%                 % 如果轮询完之后没有找到，那么就整个放弃这一帧
%                 if abs(checkCnt) == timestampShift
%                     continue
%                 end
%             end
%             disp(['lidar filename: ',lidar_dir_struct{i}]);
%             disp(['radar filename: ',radar_dir_struct{i+checkCnt}]);
        %% 读取激光雷达数据
            fid = fopen(lidar_filename);
            lidarData_frame = fread(fid);
            fclose(fid);
            lidarData_mat = zeros(606,300);
            for j=1:300
                lidarData_mat(:,j) = lidarData_frame(1+606*(j-1):606*j);
            end
            lidarData_mat = lidarData_mat(3:end,:);
            lidarDataFrame = squeeze(packLidarDataIntoFrames_v03(lidarData_mat')); % 角度*距离*帧号
%             [ranges,angles] = size(lidarDataFrame,2);
            if size(lidarDataFrame,2) < 100
                continue
            end
%             figure(4);imagesc(lidarDataFrame(:,:,1));
            nextCnt = 1;
        %% 读取毫米波数据
            % 配置参数
            if strcmp(modeButtonGroup.SelectedObject.String, 'urr')
                numADCSamples = param.numADCSamples;
                numChirpsPerFrame = param.numChirpsPerFrame;
                RX_num = param.RX_num;
                TX_num = param.TX_num;
                numDataPerFrame = numChirpsPerFrame * RX_num * numADCSamples * 2;
%                 disp(modeButtonGroup.SelectedObject.String);
            elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr')
                numADCSamples_vel = param.numADCSamples_vel;
                numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
                numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
                RX_num = param.RX_num;
                TX_num = param.TX_num;
                numDataPerFrame_v1 = numCirpsPerFrame_vel1 * RX_num * numADCSamples_vel * 2;
                numDataPerFrame_v2 = numCirpsPerFrame_vel2 * RX_num * numADCSamples_vel * 2;
                numDataPerFrame = numDataPerFrame_v1 + numDataPerFrame_v2;
                disp('mrr');    
            elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr+urr')
                numADCSamples = param.numADCSamples;
                numADCSamples_vel = param.numADCSamples_vel;
                numChirpsPerFrame = param.numChirpsPerFrame;
                numCirpsPerFrame_vel1 = param.numCirpsPerFrame_vel1;
                numCirpsPerFrame_vel2 = param.numCirpsPerFrame_vel2;
                RX_num = param.RX_num;
                TX_num = param.TX_num;
                numDataPerFrame_cp = numChirpsPerFrame * RX_num * numADCSamples * 2;
                numDataPerFrame_v1 = numCirpsPerFrame_vel1 * RX_num * numADCSamples_vel * 2;
                numDataPerFrame_v2 = numCirpsPerFrame_vel2 * RX_num * numADCSamples_vel * 2;
                numDataPerFrame = numDataPerFrame_cp + numDataPerFrame_v1 + numDataPerFrame_v2;
                disp('mrr+urr');
            end
            frame_num = 1;
            if strcmp(modeButtonGroup.SelectedObject.String, 'urr')
                disp('urr');
                radarDataAll=read2DCA1000(radar_filename, numADCSamples, RX_num, TX_num, frame_num);  %读取文件
                size(radarDataAll)
            elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr')
                disp('mrr');
                [vel1data, vel2data]=read4DCA1000(radar_filename, numADCSamples_vel, numCirpsPerFrame_vel1, numCirpsPerFrame_vel2, RX_num, frame_num);  %读取文件
                size(vel1data)
                size(vel2data)
            elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr+urr')
%                 disp('mrr+urr');
                [radarDataAll, vel1data, vel2data]=read5DCA1000(radar_filename, numADCSamples, numADCSamples_vel, numChirpsPerFrame, numCirpsPerFrame_vel1, numCirpsPerFrame_vel2, RX_num, TX_num, frame_num);  %读取文件
%                 [radarDataAll, vel1data, vel2data]=read3DCA1000('1630927957.415791.bin', 512,256,96, 128, 128, 4,3, 1);  %读取文件
%                 size(radarDataAll)
%                 size(vel1data)
%                 size(vel2data)
            end
            %% 读取视频
            frame_index = floor((i - 1) / 2) * 5 + 2 + rem(i, 2);
            frame = read(vid, frame_index);
            figure(2);subplot(2,2,2);imshow(frame);
            %% 数据处理
            linNum=1;
            if linNum > 1
                lines = 3;
            else
                lines = 1;
            end
            for lineId = 1:lines
                lidarDataFrame_singL = squeeze(lidarDataFrame(2:end,:,lineId,:));
%                 fprintf("这是第%d条线\n",lineId);

                switch modeButtonGroup.SelectedObject.String

                    case 'urr'
                        switch methodButtonGroup1.SelectedObject.String

                            case 'normal'
                                methodSign = 1;
                                [ distanceCoor, velocityCoor, distance, velocity, FinalResult, mmwavedata,dopplerSum, pcStrc ]=mmwaveResults_v1_3(radarDataAll, lidarDataFrame_singL, num2str(i), methodSign, lineId);
                            case 'doppler comp'
                                methodSign = 2;
                            case 'tracking'
                                methodSign = 3;

                        end

                    case 'mrr'
                        switch methodButtonGroup1.SelectedObject.String

                            case 'normal'
                                methodSign = 1;
                                [ distanceCoor, velocityCoor, distance, velocity, mmwavedata, dopplerSum ] = mmwaveResults_CRT_mrr(vel1data, vel2data, num2str(i));
                            case 'doppler comp'
                                methodSign = 2;
                            case 'tracking'
                                methodSign = 3;

                        end

                    case 'mrr+urr'
                        switch methodButtonGroup2.SelectedObject.String

                            case 'normal'
                                methodSign = 4;
                                [ distanceCoor_vel, velocityCoor1,velocityCoor2,distanceCoor,velocityCoor, distance, velocity, CFAROut, mmwavedata,dopplerSum, dopplerSum1,dopplerSum2,pcStrc ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(1),lineId);
                            case 'CRT'
                                methodSign = 5;
                                [ distanceCoor, velocityCoor, distance, velocity, mmwavedata, dopplerSum ] = mmwaveResults_CRT(radarDataAll, vel1data, vel2data, num2str(i));
%                                 [~] = dataCubeProcess(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(i),lineId);
                            case 'apply'
%                                 [ distanceCoor_vel, velocityCoor1,velocityCoor2,distanceCoor,velocityCoor, distance, velocity, CFAROut, mmwavedata,dopplerSum, dopplerSum1,dopplerSum2,pcStrc ] = mmwaveResults_urrsrrNormal_applied(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(1),lineId);
                                mmwaveResults_urrsrrNormal_applied(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(1),lineId);
%                                 startNum = 40;
%                                 lidarAngleGrid = (lidarData_frame(1:end-4,1)*256+lidarData_frame(1:end-4,2))/100-90;
%                                 lidarData = lidarData_frame(1:end-4,startNum:end);
%                                 [m,n] = size(lidarData);
%                                 lidarRangeGrid = 0.15*([1:n]);
%                                 figure(3);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData');
%                                 xlabel('角度');
%                                 ylabel('m');
%                                 set(gca,'YDIR','normal');
%                                 pause(0.05);
                        end

                end
                
            end
        end
        
    end
    % 删除当前包，读下一包到内存里
    function frame_change_next_Callback(~,~)
        clear radarDataAll vel1data vel2data
        nextCnt = nextCnt+1;
        fprintf("当前读取第%d包\n",nextCnt);
        len = length(filename_pack);
%         len = 1;
        if nextCnt < len
            nexFilename = filename_pack{nextCnt};
            frameNum_thisPack = param.packLength;
            [radarDataAll, vel1data, vel2data]=read3DCA1000(nexFilename, param.numADCSamples, param.numADCSamples_vel, param.numChirpsPerFrame, param.numCirpsPerFrame_vel1, param.numCirpsPerFrame_vel2, param.RX_num, param.TX_num, frameNum_thisPack);  %读取文件
        elseif nextCnt == len
            nexFilename = filename_pack{nextCnt};
            frameNum_thisPack = frame_num-(nextCnt-1)*param.packLength;
%             nexFilename = 'data.bin';
%             frameNum_thisPack = 1;
            [radarDataAll, vel1data, vel2data]=read3DCA1000(nexFilename, param.numADCSamples, param.numADCSamples_vel, param.numChirpsPerFrame, param.numCirpsPerFrame_vel1, param.numCirpsPerFrame_vel2, param.RX_num, param.TX_num, frameNum_thisPack);  %读取文件
        else
            fprintf("没有更多包了\n");
            return;
        end
        fprintf("读取结束\n");
    end

    function frame_change_stop_Callback(~,~)
        stopFlag = true;
    end
    function frame_change_reset_Callback(~,~)
        nextCnt = 0;
        fprintf("重新读第一包\n");
        frame_change_next_Callback();
    end
end