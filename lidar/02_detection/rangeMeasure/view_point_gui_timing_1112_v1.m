function [t_data,area_data,wide_data]=view_point_gui_timing_1112_v1(data,zeropoint,area_cor_form,x_cfd_cor_form,t_offset)
% ��һ֡�����ݶ�Ӧ��һ֡�Ƕ�

% load(['F:\����\',filename,'.mat']);
% a=double(a);

s=size(data,2);
t_data=zeros(s,1);
wide_data=zeros(s,1);
area_data=zeros(s,1);

% start_x=60;

for j=1:s
    %% ԭ���Ĳ�����Ѱ��ÿ�����������ֵ����������������
    % ���ڸ�һ�£����CFAR��ȡ���ķ�ֵ�����������������������Ҫ���������ʣ�����Ϊ���ģ�
    [max_data,max_index]=max(data(:,j));
%     biData = cfar_ca1D_square_dythresh(data(:,j),30,15,0.2,0);
    
    % %     t_data(j)=valve_ad_828(data(:,j),70)-58.96-0.386;
    %     %     t_data(j)=valve_ad(data(:,j),20)-55;
    %     if max_index-start_x<=3         %%�����ֵ��С���ķ�ֵ������������һ�ػز�
    %         for i=max_index:length(data(:,j))
    %             if data(i,j)<20
    %                 break
    %             end
    %         end
    %         [max_data,max_index]=max(data(i+2:end,j));  %% i+2Ϊ����һ�ػز�����ʼλ��
    %         max_index=max_index+i+2-1;
    %     end
    if max_data>5 && ~isnan(t_data(j))
        try
            [area_data(j),wide_data(j)]=area_amend_904(data(:,j),max_index);
        catch
            continue
        end
    else
        area_data(j)=NaN;
        wide_data(j)=NaN;
    end
    
    t_data(j)=cfd_ad_ys_828(data(:,j),max_index,0.6,8)+t_offset;
    
    t_cor=interp1(area_cor_form,x_cfd_cor_form,area_data(j),'linear','extrap');
%     t_cor=0;
    
    t_data(j)=t_data(j)-t_cor;
%     [~,index_zp]=max(zeropoint);
%     [area_zp,wide_zp]=area_amend_904(zeropoint,index_zp);
%     t_zp=cfd_ad_ys_828(zeropoint,index_zp,0.6,8);
%     if ~isnan(area_zp)
%         zp_cor=interp1(area_cor_form,x_cfd_cor_form,area_zp,'linear','extrap');
%     else
%         zp_cor=0;
%     end
%     
% %     zp_cor=0;
%     t_zp=t_zp-zp_cor;
%     t_data(j)=t_data(j)-t_zp;
    t_data(j)=t_data(j)-63;
end
end