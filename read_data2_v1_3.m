%%% Created by 李嘉宝 2021.05.06 version 1.2 %%%

function read_data2_v1_3
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
                  'String','normal',...
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
f.Name = '数据显示';
movegui(f,'center')
f.Visible = 'on';
global dir_path
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

    function file_button_Callback(~,~)
        [filename,dir_path] = uigetfile('.mat','打开文件');
        fprintf("开始读激光雷达数据......\n");
        filename
        data = load([dir_path,filename]);
        lidarData = data.data;
%         lidarDataFrame = packLidarDataIntoFrames(lidarData); % 角度*距离*帧号
        lidarDataFrame = packLidarDataIntoFrames_v02(lidarData); % 角度*距离*帧号
        size(lidarDataFrame)
        fprintf("结束!\n");
        fprintf("读取修正参数......\n");
        %% 载入修正参数
        area_t = load('./lidar/02_detection/rangeMeasure/measureParam/4_area_t.mat');%载入通道4修正系数
        area_cor_form4=area_t.area_t(1,:);
        x_cfd_cor_form4=area_t.area_t(2,:);
        fprintf("结束!\n");
        % load('.\measureParam\5_area_t.mat');%载入通道5修正系数
        % area_cor_form5=area_t(1,:);
        % x_cfd_cor_form5=area_t(2,:);
        % load('.\measureParam\6_area_t.mat');%载入通道6修正系数
        % area_cor_form6=area_t(1,:);
        % x_cfd_cor_form6=area_t(2,:);

    end

     function mmwavefile_Callback(~,~)
        % 配置参数
        if strcmp(modeButtonGroup.SelectedObject.String, 'urr')
            numADCSamples = param.numADCSamples;
            numChirpsPerFrame = param.numChirpsPerFrame;
            RX_num = param.RX_num;
            TX_num = param.TX_num;
            numDataPerFrame = numChirpsPerFrame * RX_num * numADCSamples * 2;
            disp(modeButtonGroup.SelectedObject.String);
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
        
        [mmwavefile_name,dir_path]=uigetfile('.bin','打开文件');        
        fprintf("开始读毫米波雷达数据\n");
        nextCnt=1;
        fprintf("当前读取第%d包\n",nextCnt);
        % read .bin file
        fid = fopen([dir_path,mmwavefile_name],'r')
        fileData = fread(fid, 'int16');
        fileSize = size(fileData, 1);
        fileBase = 1;
        adcBase = 1;
        adcData = zeros(numDataPerFrame, 1);

        while fileBase <= fileSize

            fileLength = fileData(fileBase + 2) / 2;
            adcData(adcBase: (adcBase + fileLength - 1)) = fileData((fileBase + 7): (fileBase + 7 + fileLength - 1));
            adcBase = adcBase + fileLength;
            fileBase = fileBase + 7 + fileLength;

        end

        frame_num = floor(size(adcData, 1) / numDataPerFrame)
        adcData_tail = adcData((frame_num * numDataPerFrame + 1): end);
        %tail_chirps_num = floor(length(adcData_tail) / 2 / numADCSamples);
        adcData = adcData(1: (frame_num * numDataPerFrame));
        clength = param.packLength;
        if frame_num > clength
%             res = mod(frame_num,clength);
            div = floor(frame_num/clength)
            for i=1:div
                filename = ['data_',num2str(i),'.bin'];
                filename_pack{end+1} = filename;
                adcData_temp = adcData(((i-1)*clength * numDataPerFrame)+1: (i*clength * numDataPerFrame));
                write = fopen(filename, 'wb');
                fwrite(write, adcData_temp, 'int16');
                fclose(write);
                clear adcData_temp
            end
            filename = ['data_',num2str(i+1),'.bin'];
            filename_pack{end+1} = filename;
            adcData_temp = adcData((div*clength * numDataPerFrame)+1: (i*frame_num * numDataPerFrame));
            write = fopen(filename, 'wb');
            fwrite(write, adcData_temp, 'int16');
            fclose(write);
            clear adcData_temp
        else
            filename = 'data.bin';
            write = fopen(filename, 'wb');
            fwrite(write, adcData, 'int16');
            fclose(write);
        end

        writetail = fopen('data_tail.bin', 'wb');
        fwrite(writetail, adcData_tail, 'int16');
        fclose(writetail);
        clear adcData
        clear adcData_tail
        sta = fclose(fid);
        
        if strcmp(modeButtonGroup.SelectedObject.String, 'urr')
            disp('urr');
            radarDataAll=read2DCA1000('data.bin', numADCSamples, RX_num, TX_num, frame_num);  %读取文件
            fidSize = size(radarDataAll)
        elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr')
            disp('mrr');
            [vel1data, vel2data]=read4DCA1000('data.bin', numADCSamples_vel, numCirpsPerFrame_vel1, numCirpsPerFrame_vel2, RX_num, frame_num);  %读取文件
            fidSize1 = size(vel1data)
            fidSize2 = size(vel2data)
        elseif strcmp(modeButtonGroup.SelectedObject.String, 'mrr+urr')
            disp('mrr+urr');
            if length(filename_pack) > 1
                disp('分包读第一包');
                filename = filename_pack{1};
                [radarDataAll, vel1data, vel2data]=read3DCA1000(filename, numADCSamples, numADCSamples_vel, numChirpsPerFrame, numCirpsPerFrame_vel1, numCirpsPerFrame_vel2, RX_num, TX_num, clength);  %读取文件
            else
                [radarDataAll, vel1data, vel2data]=read3DCA1000('data.bin', numADCSamples, numADCSamples_vel, numChirpsPerFrame, numCirpsPerFrame_vel1, numCirpsPerFrame_vel2, RX_num, TX_num, frame_num);  %读取文件
%                 [radarDataAll, vel1data, vel2data]=read3DCA1000('1630927957.415791.bin', 512,256,96, 128, 128, 4,3, 1);  %读取文件
            end
            fidSize = size(radarDataAll)
            fidSize1 = size(vel1data)
            fidSize2 = size(vel2data)
        end
%         global vel1data
%         global vel2data
%         fid2=readTailDCA1000('data_tail.bin');

        fprintf("结束!\n");
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
%         [angle,range,frame] = size(lidarDataFrame);
        [angle,range,linNum,frame] = size(lidarDataFrame);        
        stopFlag = false;
        frameIdx = str2num(frameIndex.String);
        Temp1_accumulate = [];
%         figure(2);
        for i=frameIdx:frame
            if stopFlag==true
                fprintf("停止\n");
                break;
            end
            this_frameNum = i+(nextCnt-1)*param.packLength;
            fprintf("这是第%d帧\n",this_frameNum);
%             pcStrcAll = [];
            all_x = [];
            all_y = [];
            all_z = [];
            all_i = [];
            linNum=1;
            if linNum > 1
                lines = 3;
            else
                lines = 1;
            end
            stopFlag = false;
            for lineId = 1:lines
                if stopFlag
                    fprintf('停止!\n');
                    break
                end
                lidarDataFrame_singL = squeeze(lidarDataFrame(:,:,lineId,:));
                lidarData_frame = lidarDataFrame_singL(:,:,this_frameNum); % 取第i帧激光雷达
                fprintf("这是第%d条线\n",lineId);
%                 figure(3);imagesc(lidarData_frame(1:end-4,14:end)');
%                 set(gca,'YDIR','normal');
%                 pause(0.01);
%               [~]=mmwaveResults_v1_2(radarDataAll, lidarDataFrame, num2str(i));

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
%                                 [ distanceCoor_vel, velocityCoor1,velocityCoor2,distanceCoor,velocityCoor, distance, velocity, CFAROut, mmwavedata,dopplerSum, dopplerSum1,dopplerSum2,pcStrc ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(i),lineId,Temp1_accumulate);
                                [ distanceCoor_vel, velocityCoor1,velocityCoor2,distanceCoor,velocityCoor, distance, velocity, CFAROut, mmwavedata,dopplerSum, dopplerSum1,dopplerSum2,pcStrc ] = mmwaveResults_urrsrrNormal(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(i),lineId,Temp1_accumulate);
                            case 'CRT'
                                methodSign = 5;
                                [ distanceCoor, velocityCoor, distance, velocity, mmwavedata, dopplerSum ] = mmwaveResults_CRT(radarDataAll, vel1data, vel2data, num2str(i));
%                                 [~] = dataCubeProcess(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(i),lineId);
                            case 'apply'
                                mmwaveResults_urrsrrNormal_applied(radarDataAll, vel1data, vel2data,lidarDataFrame_singL, num2str(i),lineId,Temp1_accumulate);
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