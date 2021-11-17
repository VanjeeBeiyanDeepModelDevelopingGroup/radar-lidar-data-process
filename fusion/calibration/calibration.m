% function [ret] = calibration(dirPath)
%CALIBRATION 联合标定脚本，生成标定外参数
%   点对点之后再进行最小二乘优化即可得到标定参数
dirPath = './fusion/calibration/calibData/';
dirContent = dir(dirPath);
num = length(dirContent);
oldtxt = 28;
figure(8);
for i = 3:num
    filename = dirContent(i).name;
    filename = [dirPath,filename];
    result = load(filename);
    % 读到的数据是radar angle, radar range, lidar angle, lidar range
    result = result.result;
    color = [rand,rand,rand];
    plot(result(1),result(2),'.','color',color,'markersize',40);hold on
    plot(result(3),result(4),'p','color',color,'markersize',40);
    line([result(1),result(3)],[result(2),result(4)],'color',color);
    txt = filename(end-21:end-20);
    text(result(1),result(2),num2str(txt),'fontsize',35);
    xlim([-20,20]);
    ylim([10,15]);
    xlabel('角度');
    ylabel('距离');
end
legend('毫米波','激光');



% end
%% 对比一下测量角度与真实角度
true_angles = [-45:5:-10,0:5:45];
test_angles = [-52,-53.5,-39.333,-41,-5.75,-8.25,-9.5,-10,2,7,10.5,11,11,19.75,19.67,20,59.5,60.5];
figure;subplot(2,1,1);
plot(true_angles,'ro-');hold on 
plot(test_angles,'bp-');
legend('真实值','测量值');
subplot(2,1,2);
plot(true_angles-test_angles,'p');
% 只取-10到10度的数据进行校正
calibTrue_angles = [-10,0:5:10];
calibTest_angles = [-3.833,2,7,10.5];
delta_angles = calibTrue_angles-calibTest_angles;
% 结果
true = -0.0255*test^2+1.527*test-3.596;