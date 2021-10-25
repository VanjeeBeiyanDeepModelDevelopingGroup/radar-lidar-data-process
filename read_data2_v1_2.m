%%% Created by 李嘉宝 2021.02.19 version 1.2 %%%
%%% 增加距离，速度，xyz坐标列表输出 %%%



function read_data2_v1_2
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
htext = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[180,560,200,30],...
    'Background','w','HorizontalAlignment','left');
rangetext = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[430,560,200,30],...
    'Background','w','HorizontalAlignment','left');
mmwavefile = uicontrol('Style','pushbutton','Units','pixels',...
    'String','毫米波文件','FontSize',12,...
    'Position',[50,516,100,30],...
    'Callback',@mmwavefile_Callback);
frametext = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[177,510,200,30],...
    'String','请输入想查看的帧索引：','HorizontalAlignment','left');
frameIndex = uicontrol('Style','edit','Units','pixels',...
    'String','','FontSize',12,...
    'Position',[357,515,35,30]);
frameNumtext = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[668,555,200,30],...
    'Visible','off','HorizontalAlignment','left');
framePlotButton = uicontrol('Style','pushbutton','Units','pixels',...
    'String','画图','FontSize',12,...
    'Position',[430,516,50,30],...
    'Callback',@frame_change_disp_Callback);
modeButtonGroup = uibuttongroup('Visible','on','Units','pixels',...
                  'Position',[500 516 120 40],'Title','Modes',...
                  'SelectionChangedFcn',@bselection);
mode1 = uicontrol(modeButtonGroup,'Style','radiobutton',...
                  'String','urr',...
                  'Position',[10 10 30 15],...
                  'HandleVisibility','on');
mode2 = uicontrol(modeButtonGroup,'Style','radiobutton',...
                  'String','urr+srr',...
                  'Position',[60 10 50 15],...
                  'HandleVisibility','on'); 
methodButtonGroup1 = uibuttongroup('Visible','off','Units','pixels',...
                  'Position',[650 516 105 80],'Title','Methods',...
                  'SelectionChangedFcn',@m1selection);
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
                  'Position',[650 516 105 80],'Title','Methods',...
                  'SelectionChangedFcn',@m2selection);
method1 = uicontrol(methodButtonGroup2,'Style','radiobutton',...
                  'String','normal',...
                  'Position',[10 40 50 15]);
method2 = uicontrol(methodButtonGroup2,'Style','radiobutton',...
                  'String','CRT',...
                  'Position',[10 15 50 15]);
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
finalResultTable.Units = 'normalized';
modeButtonGroup.Units = 'normalized';
methodButtonGroup.Units = 'normalized';
f.Units = 'normalized';
f.Name = '数据显示';
movegui(f,'center')
f.Visible = 'on';
global dir_path
global radarDataAll
global mmwavefile_name
% global lidarData;
global lidarDataFrame
    function file_button_Callback(~,~)
        [filename,dir_path] = uigetfile('.mat','打开文件');
        fprintf("开始读激光雷达数据\n");
        data = load([dir_path,filename]);
        lidarData = data.data;
        lidarDataFrame = packLidarDataIntoFrames(lidarData); % 角度*距离*帧号
        fprintf("结束!\n");

    end

     function mmwavefile_Callback(~,~)
        % 配置参数
        numADCSamples = param.numADCSamples;
        numChirpsPerFrame = param.numChirpsPerFrame;
        RX_num = param.RX_num;
        TX_num = param.TX_num;
        numDataPerFrame = numChirpsPerFrame * RX_num * numADCSamples * 2;
        
        mmwavefile_name=uigetfile('.bin','打开文件');        
        fprintf("开始读毫米波雷达数据\n");
        % read .bin file
        fid = fopen(mmwavefile_name,'r');
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

        frame_num = floor(size(adcData, 1) / 2 / numADCSamples / numChirpsPerFrame / RX_num)
        adcData_tail = adcData((frame_num * 2 * numADCSamples * numChirpsPerFrame + 1): end);
        tail_chirps_num = floor(length(adcData_tail) / 2 / numADCSamples);
        adcData = adcData(1: (frame_num * 2 * numADCSamples * numChirpsPerFrame * RX_num));

        write = fopen('data.bin', 'wb');
        fwrite(write, adcData, 'int16');
        fclose(write);

        writetail = fopen('data_tail.bin', 'wb');
        fwrite(write, adcData_tail, 'int16');
        fclose(write);

        sta = fclose(fid);
        % frame_num = floor(frame_num/10)
        radarDataAll=read2DCA1000('data.bin', numADCSamples, RX_num, TX_num, frame_num);  %读取文件
        % save fid fid
        fidSize = size(radarDataAll)
%         fid2=readTailDCA1000('data_tail.bin');
        fprintf("结束!\n");
     end

    function bselection(~,event)
        if strcmp(event.NewValue.String, 'urr')
            methodButtonGroup2.Visible = 'off';
            methodButtonGroup1.Visible = 'on';
        elseif strcmp(event.NewValue.String, 'urr+srr')
            methodButtonGroup1.Visible = 'off';
            methodButtonGroup2.Visible = 'on';
        end
    end

    function frame_change_disp_Callback(~,~)
        [angle,range,frame] = size(lidarDataFrame);
        % 画出激光雷达波形，找到对应帧
        for i=1:frame
%             lidarData_frame = lidarDataFrame(:,:,str2num(frameIndex.String)); % 取第i帧激光雷达
%             fprintf("这是第%d帧\n",str2num(frameIndex.String));
            lidarData_frame = lidarDataFrame(:,:,i); % 取第i帧激光雷达
            fprintf("这是第%d帧\n",i);
%             [~]=mmwaveResults_v1_2(radarDataAll, lidarDataFrame, num2str(i));

            [ distanceCoor, velocityCoor, distance, velocity, FinalResult, mmwavedata, dopplerSum ]=mmwaveResults_v1_2(radarDataAll, lidarDataFrame, num2str(i));
%             [ distanceCoor, velocityCoor, distance, velocity, FinalResult, mmwavedata, dopplerSum ]=mmwaveResults_v1_2(radarDataAll, lidarDataFrame, frameIndex.String);
            figure(2);
%             cla reset;
            clf();
            subplot(2,1,1);
%             lidarDisCor = 0:0.15:max(distanceCoor);
%             lidarIndex = 1: length(lidarDisCor);
%              for j=150:180
             for j=160:161
                 lidarDataString = lidarData_frame(j,60:603)*100;
                 lidarDataString = lidarDataString-mean(lidarDataString);
%                  plot(distanceCoor,mmwavedata(1,:,1,1),lidarDisCor,lidarDataString(lidarIndex));hold on
                lidarDisCor = 0:0.15:((size(lidarDataString, 2)-1)*0.15);
                plot(distanceCoor,mmwavedata(1,:,1,1),lidarDisCor,lidarDataString);hold on
             end
            subplot(2,1,2);imagesc(distanceCoor,velocityCoor,dopplerSum');title('单个天线的DopplerFFT波形');
            set(gca,'YDir','normal');
            xlabel('距离(m)');ylabel('速度(m/s)');% zlabel('频谱幅值');
%             finalResultTable.Data = FinalResult;
%             axes(plot_data);
% %             clf(plot_data);
%             plot(distanceCoor,mmwavedata(1,:,1,1));hold on
% %              for j=150:180
% %                  lidarDataString = lidarData_frame(j,:)*20;
% %                  plot(lidarDataString);hold on
% %              end
%             axes(plot_mmwave_data);
%             plot(distance, velocity, 'rp');
%             axis([min(distanceCoor) max(distanceCoor) min(velocityCoor) max(velocityCoor)]);
% 
%             finalResultTable.Data = FinalResult;
% 
%             pc = FinalResult(:,3:5);
%             x = pc(:,1);
%             y = pc(:,2);
%             z = pc(:,3);
        end
    end

end
