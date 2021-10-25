function [pcStrc] = calElevAngLidar(line_id,allMiuSig)
%calElevAngLidar 计算多线的垂直角度
%   1-8线基本上是水平的，9-16线是歪的
%   输入是每条线扫描后融合输出的点，输出是三维空间的点
% [angleNum,adNum] = size(allMiuSig);
angle = allMiuSig(:,1);
t_data = allMiuSig(:,2);
line_id = line_id+1;
if line_id>=1 && line_id<=8
    ver_angle= (-2.17 + 0.62*(line_id-1))*ones(size(angle));
else
    w=[-3.4248338138919e-07,-2.28509524735483e-05,0.0396487190040173,2.14160921236336];
    ang=w(1)*angle.^3 + w(2)*angle.^2 + w(3)*angle + w(4);
    ver_angle=ang-(0.62*(15-(line_id-1)));
end
B=ver_angle*pi/180;
A=angle*pi/180;
pcStrc.vertex.x=-t_data.*cos(A).*cos(B);
pcStrc.vertex.y=t_data.*sin(A).*cos(B);
pcStrc.vertex.z=t_data.*sin(B);
% pcStrc.vertex.i = allMiuSig(:,3);
% 显示
% figure(2);
% scatter3(pcStrc.vertex.x,pcStrc.vertex.y,pcStrc.vertex.z,50,'.');view([0,0,1]);
% xlim([-10,10]);ylim([0,80]);title("点云");
% grid off;axis on;
end