function view_gui_read_1112
%��һ֡���ݶ�Ӧ��һ֡�Ƕ�

close all
%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[0,30,1370,650]);

%  Construct the components.
hconf = uicontrol('Style','pushbutton','Units','pixels',...
    'String','��ʼ','FontSize',12,...
    'Position',[50,600,70,30],...
    'Callback',@conf_button_Callback);

channel=uicontrol(f,'Style','text','Units','pixels');
channel.Position=[1000,600,100,30];
channel.String='�����ٶ�';
channel.FontSize=12;
hline = uicontrol('Style','togglebutton','String','����',...
    'FontSize',12,'Position',[250,600,70,30]);
hcomp = uicontrol('Style','togglebutton','String','�ȶ�',...
    'FontSize',12,'Position',[350,600,70,30]);
hclear = uicontrol('Style','pushbutton',...
    'FontSize',12,'String','���',...
    'Position',[500,600,70,30],...
    'Callback',@clear_button_Callback);
hspeed = uicontrol('Style','edit','String','1',...
    'FontSize',12,'Position',[1100,600,30,30]);
hwait = uicontrol('Style','pushbutton',...
    'FontSize',12,'String','��ͣ',...
    'Position',[650,600,70,30],...
    'Callback',@wait_button_Callback);
hstop = uicontrol('Style','pushbutton',...
    'FontSize',12,'String','ֹͣ',...
    'Position',[800,600,70,30],'Callback',@stop_button_Callback);
hrec = uicontrol('Style','pushbutton',...
    'FontSize',12,'String','�ط�',...
    'Position',[900,600,70,30],'Callback',@rec_button_Callback);
point_cloud = axes('Units','Pixels','Position',[50,140,400,400]);
wave_plot=axes('Units','Pixels','Position',[500,140,400,400]);
wide_area_plot=axes('Units','Pixels','Position',[950,140,400,400]);
htext2 = uicontrol(f,'Style','text','Units','pixels',...
    'FontSize',12,'Position',[50,20,1000,70],...
    'Background','w','HorizontalAlignment','left');
hsave = uicontrol('Style','text',...
    'FontSize',12,'String','��1',...
    'Position',[950,600,70,30]);
hoverl = uicontrol('Style','togglebutton','String','����',...
    'FontSize',12,'Position',[1150,600,70,30]);

hslider = uicontrol('Style','slider',...
    'Position',[40,100,410,15]);
stext = uicontrol('Style','edit','Units','pixels',...
    'FontSize',12,'Position',[460,100,50,15],...
    'String','1');

align([hconf,hline,hcomp,hclear,hwait,hstop,hrec,hsave,channel,hspeed,hoverl],'Distribute','Center');
f.Units = 'normalized';
point_cloud.Units = 'normalized';
wide_area_plot.Units = 'normalized';
wave_plot.Units = 'normalized';
hconf.Units = 'normalized';
htext2.Units = 'normalized';
hline.Units = 'normalized';
hcomp.Units = 'normalized';
hwait.Units = 'normalized';
hstop.Units = 'normalized';
hrec.Units = 'normalized';
hclear.Units = 'normalized';
hspeed.Units = 'normalized';
hoverl.Units = 'normalized';
hconf.Units = 'normalized';
channel.Units = 'normalized';
hsave.Units = 'normalized';
hslider.Units = 'normalized';
stext.Units = 'normalized';
f.Name = '������ʾ';
% movegui(f,'center')
f.Visible = 'on';
wait_flag=0;
stop_flag=0;
target_x=NaN;
target_y=NaN;
fre=0;
colors=[[0,1,0]
    [1,0.5,0]
    [0.5,0.5,0.5]
    [0.5,0.2,0.2]
    [0,0,1]
    [1,0,0]
    [0,0.5,1]
    [0,1,1]
    [0,0.8,0.3]
    [1,0,1]
    [0.5,0,1]
    [0.8,0,0.3]
    [0,0,0]
    [0.2,0.4,0.6]
    [0.6,0.4,0.2]
    [0.2,0.6,0.4]];
for line_n=1:16
    line_text= uicontrol(f,'Style','text','Units','pixels');
    line_text.Position = [60*line_n-50,550,30,30];
    line_text.String=num2str(line_n);
    line_text.FontSize=12;
    line_text.Units = 'normalized';
    line_bg= uicontrol(f,'Style','text','Units','pixels');
    line_bg.Position = [60*line_n-20,565,20,20];
    line_bg.BackgroundColor = colors(line_n,:);
    line_bg.Units = 'normalized';
end
area_t = [];
% load('.\��·����\1_area_t.mat');%����ͨ��1����ϵ��
% area_cor_form1=area_t(1,:);
% x_cfd_cor_form1=area_t(2,:);
% load('.\��·����\2_area_t.mat');%����ͨ��2����ϵ��
% area_cor_form2=area_t(1,:);
% x_cfd_cor_form2=area_t(2,:);
% load('.\��·����\3_area_t.mat');%����ͨ��3����ϵ��
% area_cor_form3=area_t(1,:);
% x_cfd_cor_form3=area_t(2,:);
load('.\measureParam\4_area_t.mat');%����ͨ��4����ϵ��
area_cor_form4=area_t(1,:);
x_cfd_cor_form4=area_t(2,:);
load('.\measureParam\5_area_t.mat');%����ͨ��5����ϵ��
area_cor_form5=area_t(1,:);
x_cfd_cor_form5=area_t(2,:);
load('.\measureParam\6_area_t.mat');%����ͨ��6����ϵ��
area_cor_form6=area_t(1,:);
x_cfd_cor_form6=area_t(2,:);
% load('.\��·����\7_area_t.mat');%����ͨ��7����ϵ��
% area_cor_form7=area_t(1,:);
% x_cfd_cor_form7=area_t(2,:);
% load('.\��·����\8_area_t.mat');%����ͨ��8����ϵ��
% area_cor_form8=area_t(1,:);
% x_cfd_cor_form8=area_t(2,:);

load('.\��·����\std_wave.mat','std_wave');
line_id=0;
t_offset=0;
area_cor_form=[];
x_cfd_cor_form=[];
ui_dir=' ';
%    view_point_gui('8.4����\8')
    function conf_button_Callback(~,~)
        %         cd(['F:\����\',hfn.String])
        ui_dir=uigetdir('F:\����\11.16����');
%         for mm=1:16
        for mm=4:6
            max_c=length(dir([ui_dir,'\��',num2str(mm),'\','*.mat']))-length(dir([ui_dir,'\','*.mat']));
%             max_c=length(dir([ui_dir,'\','*.mat']));
            hslider.Min=1;
            hslider.Max=max_c;
            stext.String='1';
            target_x=NaN;
            target_y=NaN;

            line_id=mm;

            switch line_id
                case {1,9}
                    t_offset=0.4246;
                    area_cor_form=area_cor_form1;
                    x_cfd_cor_form=x_cfd_cor_form1;
                case {2,10}
                    t_offset=0.3144;
                    area_cor_form=area_cor_form2;
                    x_cfd_cor_form=x_cfd_cor_form2;
                case {3,11}
                    t_offset=0.5254;
                    area_cor_form=area_cor_form3;
                    x_cfd_cor_form=x_cfd_cor_form3;
                case {4,12}
                    t_offset=0.4893;
                    area_cor_form=area_cor_form4;
                    x_cfd_cor_form=x_cfd_cor_form4;
                case {5,13}
                    t_offset=0.4546;
                    area_cor_form=area_cor_form5;
                    x_cfd_cor_form=x_cfd_cor_form5;
                case {6,14}
                    t_offset=0.5160;
                    area_cor_form=area_cor_form6;
                    x_cfd_cor_form=x_cfd_cor_form6;
                case {7,15}
                    t_offset=0.4740;
                    area_cor_form=area_cor_form7;
                    x_cfd_cor_form=x_cfd_cor_form7;
                case {8,16}
                    t_offset=0.4973;
                    area_cor_form=area_cor_form8;
                    x_cfd_cor_form=x_cfd_cor_form8;
                otherwise
                    error('ͨ��ѡ�����')
            end
            plot_fig;
        end
    end

    function plot_fig
        
        hslider.Value=1;
        wait_flag=0;
        hwait.String='��ͣ';
        stop_flag=0;
        while 1
            if stop_flag==1
                disp('��ֹͣ')
                break
            end
            %             try
            file_num=round(hslider.Value);
            filename=[ui_dir,'\��',num2str(line_id),'\',num2str(file_num)];
            if exist([filename,'.mat'],'file')
                load([filename,'.mat'],'a')      %%��ȡ�ļ�
                data=double(a.data-126);
                angle=double((a.angle_o(1,:)*256+a.angle_o(2,:)))*0.01;
                angle=angle';
                zeropoint=double(a.zeropoint-126);
                
                if line_id==2 || line_id==4 || line_id==6 || line_id==10 || line_id==12 || line_id==14
                    for i=1:size(data,2)    %�˵�С��
                        y=data(:,i);
                        for j=2:length(y)-1
                            if (y(j)-y(j-1)>1 && y(j)-y(j+1)>1) || (y(j-1)-y(j)>1 && y(j+1)-y(j)>1)
                                y(j)=(y(j-1)+y(j+1))/2;
                            end
                        end
                        data(:,i)=y;
                    end
                    y=zeropoint;
                    for j=39:8:95
                        if (y(j)-y(j-1)>=1 && y(j)-y(j+1)>=1) || (y(j-1)-y(j)>=1 && y(j+1)-y(j)>=1)
                            y(j)=(y(j-1)+y(j+1))/2;
                        end
                    end
                    zeropoint=y;
                end
            
                [t_data,area_data,wide_data]=view_point_gui_timing_1112(data,zeropoint,area_cor_form,x_cfd_cor_form,t_offset);
                
            else
                disp('�Ѷ�ȡ���')
                break
            end
            %             catch ME
            %                 disp(ME)
            %                 break
            %             end
            if ~hoverl.Value
                cla(point_cloud)
            end
            cla(wide_area_plot)
            cla(wave_plot)
            htext2.String='';
            %         filename=hfn.String;
            
            if line_id>=1 && line_id<=8                                               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
                ver_angle= (-2.17 + 0.62*(line_id-1))*ones(size(angle));
            else
                w=[-3.4248338138919e-07,-2.28509524735483e-05,0.0396487190040173,2.14160921236336];
                ang=w(1)*angle.^3 + w(2)*angle.^2 + w(3)*angle + w(4);
                ver_angle=ang-(0.62*(15-(line_id-1)));
            end
            B=ver_angle*pi/180;
            A=angle*pi/180;
            x_point=-t_data*0.1498.*cos(A).*cos(B);
            y_point=t_data*0.1498.*sin(A).*cos(B);
            z_point=t_data*0.1498.*sin(B);
            
            scatter3(point_cloud,x_point,y_point,z_point,4,colors(line_id,:),'filled')
            %             drawnow
            title(point_cloud,'����','FontWeight','bold','FontSize',15)
            grid(point_cloud,'on')
            hold(point_cloud,'on')
            %             axis(point_cloud,'manual')
            if fre==0
                xlim(point_cloud,[-10,10])
                ylim(point_cloud,[0,20])
                zlim(point_cloud,[-10,10])
                set(point_cloud,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
                view(point_cloud,0,90)
                %                 axis(point_cloud,'equal')
            end
            axis(point_cloud,'manual')
            
            
            plot(wide_area_plot,zeropoint,'linewidth',1.5)
            hold(wide_area_plot,'on')
            if fre==0
                xlim(wide_area_plot,[0,150])
                ylim(wide_area_plot,[-30,130])
            end
            axis(wide_area_plot,'manual')
            if ~isnan(target_x)
                min_dis=inf;
                for i=1:length(x_point)
                    if (x_point(i)-target_x)^2+(y_point(i)-target_y)^2<min_dis
                        min_dis=(x_point(i)-target_x)^2+(y_point(i)-target_y)^2;
                        index=i;
                    end
                end
                if min_dis<0.001
                    scatter(point_cloud,x_point(index),y_point(index),10,'r')
                    htext2.String ={'������Ϣ:',['t:',num2str(t_data(index)),...
                        '     x:',num2str(x_point(index)),'      y:',num2str(y_point(index))]};
                    plot(wave_plot,data(:,index),'linewidth',1.5)
                end
            end
            %         axis equal
            % set(gca,'LooseInset',get(gca,'TightInset'))
            %         datacursormode off
            %             datacursormode on
            dcm_obj = datacursormode(f);
            set(dcm_obj,'UpdateFcn',{@myupdatefcn,data,t_data,angle,area_data,...
                wide_data,std_wave})
            %             file_num=file_num+1;
            %             pause(0.05)
            %             set(f,'DoubleBuffer','on')
            %             drawnow
            a=str2double(hspeed.String);
            pause(0.5/a)
            %             file_num=file_num+1;
            hslider.Value=hslider.Value+1;
            stext.String=num2str(round(hslider.Value));
            
            fre=1;
        end
    end
        
    function txt = myupdatefcn(source,event_obj,data,t_data,angle,area_data,...
            wide_data,std_wave)
        % Customizes text of data tips
        %         ind=get(event_obj,'DataIndex');
        
        h = get(point_cloud, 'children');
        if length(h)>=2
            delete(h(1));
        end
        
        pos=get(event_obj,'Position');
        index=get(event_obj,'DataIndex');
        if source.Parent==point_cloud||source.Parent==wide_area_plot
            if hline.Value==0
                cla(wave_plot)
            end
            plot(wave_plot,data(:,index),'linewidth',1.5)
            hold(wave_plot,'on')
            title(wave_plot,'����','FontWeight','bold','FontSize',15)
            if hcomp.Value==1
                cla(wave_plot)
                for i=1:length(std_wave.pw)
                    if wide_data(index)>=std_wave.pw(i)
                        data_comp=std_wave.data(:,i);
                        data_comp2=std_wave.data(:,i-1);
                        break
                    end
                end
                plot(wave_plot,(1:length(data_comp))+floor(t_data(index)-35),data_comp,'r','linewidth',1.5)
                plot(wave_plot,(1:length(data_comp))+floor(t_data(index)-35),data_comp2,'k','linewidth',1.5)
                hold(wave_plot,'on')
                plot(wave_plot,data(:,index),'b','linewidth',1.5)
            end
            xlim(wave_plot,[t_data(index)+40,t_data(index)+190])
            ylim(wave_plot,[-30,130])
            txt={['x:',num2str(pos(1))],['y:',num2str(pos(2))]};
            
            htext2.String ={'������Ϣ:',['    ���:',num2str(index),'              ʱ��: ',num2str(t_data(index)),'          �Ƕ�: ',num2str(angle(index)),...
                '              ���: ',num2str(area_data(index)),'    ����:',num2str(wide_data(index))]};
        else
            txt={['x:',num2str(pos(1))],['y:',num2str(pos(2))]};
        end
        target_x=pos(1);
        target_y=pos(2);
    end
    function clear_button_Callback(~,~)
        cla(wave_plot)
        htext2.String=' ';
        htext2.BackgroundColor=[1 1 1];
        target_x=NaN;
        target_y=NaN;
        h = get(point_cloud, 'children');
        if length(h)>=2
            delete(h(1));
        end
    end
    function wait_button_Callback(~,~)
        if wait_flag==0
            wait_flag=1;
            hwait.String='����';
            uiwait(f)
        else
            wait_flag=0;
            hwait.String='��ͣ';
            uiresume(f)
        end
    end
    function stop_button_Callback(~,~)
        stop_flag=1;
        uiresume(f)
        hwait.String='��ͣ';
    end
    function rec_button_Callback(~,~)
        plot_fig;
    end
end
