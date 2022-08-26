%% 静态误差曲线
gt_range = 5;
radar_range = 4.92;
gt_angles = [-45,-30:5:30];
% fft_angles = [-29.26,-30.73,-31.48,-32.23,-81.45,-90.0,-81.45,81.45,77.89,75.16	30.73,29.26,30.0,29.26];
music_angles_dataformat_1_4rx = [-42,-34,-23,-19,-10,-5,0,0,5,12,15,23,27,32];
music_angles_dataformat_2_4rx = [-44,-32,-22,-21,-10,-5,0,0,6,12,17,23,28,33];

music_angles_dataformat_1_8rx_13combination = [-32,-29,-28,-27,-2,-1,0,0,2,3,26,29,30,31];
music_angles_dataformat_2_8rx_13combination = [-31,-37,-20,-23,-2,-5,4,0,4,7,21,25,28,38];

music_angles_dataformat_1_8rx_31combination = [-32,-29,-27,-27,-2,-1,0,0,2,3,27,29,30,31];
music_angles_dataformat_2_8rx_31combination = [-34,-21,-13,-31,-2,4,-3,0,-1,21,10,32,31,24];

music_angles_dataformat_1_12rx = [-30,-29,-29,-28,-1,0,0,0,1,2,28,30,30,31];
music_angles_dataformat_2_12rx = [-30,-34,-23,-26,-1,-3,3,0,2,4,25,28,29,35];


% 画误差曲线
figure(4);
% subplot(2,1,1);
plot(gt_angles,gt_angles,'ro-');hold on
plot(gt_angles,music_angles_dataformat_1_4rx,'o-');hold on
plot(gt_angles,music_angles_dataformat_2_4rx,'o-');hold on
plot(gt_angles,music_angles_dataformat_1_8rx_13combination,'o-');hold on
plot(gt_angles,music_angles_dataformat_2_8rx_13combination,'o-');hold on
plot(gt_angles,music_angles_dataformat_1_8rx_31combination,'o-');hold on
plot(gt_angles,music_angles_dataformat_2_8rx_31combination,'o-');hold on
plot(gt_angles,music_angles_dataformat_1_12rx,'o-');hold on
plot(gt_angles,music_angles_dataformat_2_12rx,'o-');hold on
% legend('gt','fft','music');
legend('gt','间隔读-4RX','顺序读-4RX','间隔读13组合','顺序读13组合','间隔读31组合','顺序读31组合','间隔读-12RX','顺序读-12RX');
% 画点分布位置图
figure(5);
for i = 1:length(gt_angles)
    gt_angle = gt_angles(i);
    gt_pos = [gt_range*sin(gt_angle/180*pi), gt_range*cos(gt_angle/180*pi)];
    plot(gt_pos(1),gt_pos(2),'ro');hold on
%     fft_angle = fft_angles(i);
%     fft_pos = [radar_range*sin(fft_angle/180*pi), radar_range*cos(fft_angle/180*pi)];
%     plot(fft_pos(1), fft_pos(2),'go');hold on
    music_angle = music_angles(i);
    music_pos = [radar_range*sin(music_angle/180*pi), radar_range*cos(music_angle/180*pi)];
    plot(music_pos(1), music_pos(2),'bo');hold on
end
% legend('ground truth','fft','music');
legend('ground truth','music');
hold off;