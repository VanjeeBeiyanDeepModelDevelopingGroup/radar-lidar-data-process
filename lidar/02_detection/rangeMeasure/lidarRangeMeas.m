function [pcStrc,allMiuSig] = lidarRangeMeas(lidarData,line_id)
% LIDARRANGEMEAS ׼ȷ��࣬��������
% lidarData: ���ߵĸ��Ƕ�AD����
%% ����������ȫ�ֱ���
global area_cor_form4
global x_cfd_cor_form4
%% �������ݴ���
angle = (lidarData(1,:)*256+lidarData(2,:))*0.01;
lidarData = lidarData(4:end,:)-127;
%% ��������
[ads,angles] = size(lidarData); % ADnum*angles
<<<<<<< HEAD
% t_offset = 0.4546;
t_offset = 11.5226;
=======
t_offset = 11.5226;
% t_offset = 0;
>>>>>>> 0ddd06f38478382d6bb7885d3553995dd4b8069f
zeropoint = 2*ones(1,angles);
[t_data,area_data,wide_data]=view_point_gui_timing_1112_v1(lidarData,zeropoint,area_cor_form4,x_cfd_cor_form4,t_offset);

%% �����������xyz
allMiuSig = [];
allMiuSig(:,1) = angle;
allMiuSig(:,2) = 0.15*t_data;
pcStrc = calElevAngLidar(line_id,allMiuSig);

end

