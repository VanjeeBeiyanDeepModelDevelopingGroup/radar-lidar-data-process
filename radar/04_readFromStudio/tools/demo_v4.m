%
% 简介：
% 利用ti给的rawDataReader.m与verify_data.m
% 来进行studio数据的读取
% 注意：需要进行TDM三发四收设置
% 
% 操作步骤：
% 1. 先利用mmwave studio进行数据的采集，具体来说，需要设置studio进行TDM三发四收
% 1.1 profile不用变，TX1,TX2,TX3使能，chirp 00 11 22，frame 0-2
% 1.2 具体图见这个链接：https://e2e.ti.com/support/sensors-group/sensors/f/sensors-forum/713215/dca1000evm-awr1642evm-mmwave-studio-tdm-or-bpm-mimo-setting
%
% 2. 保存raw_adc.bin，xx.setup.json，mmwave.json文件。其中前者为数据文件，后二者为当前数据的配置文件
% 2.1 注意修改setup.json中的configUsed和fileBasePath到数据对应目录下，数据名称processedFileName修改保存名称
%
% 3. 将数据与配置文件放入data文件夹，进行数据可视化
%
%
% update: 2022-09-08，demo_v4实现双子帧功能
%
%
close all
clear all
c = 3e8; % m/s
% setup_filename = './data/0816/1.setup.json';
setup_filename = './data/0905_TI外场/配置/MIMO_TX0&TX2/1.setup.json';
%% 可调参数
range_cfar_thresh = 0.25;
doppler_cfar_thresh = 0.1;
frame_number = 1;
custome_peak_indx = 1;  % 调整这个值选取第几个峰为目标峰
MIMO = 0;
%% 开始
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
frame = radar_cube_struct.data{frame_number};
if(MIMO)
    radar_cube_tx1 = frame(1:3:end,:,:);
    radar_cube_tx2 = frame(2:3:end,:,:);
    radar_cube_tx3 = frame(3:3:end,:,:);
    radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);
else
    radar_cube_all = frame;
end

[chirps,atennas,range_all_indx] = size(radar_cube_all);
range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(0:range_all_indx-1);
vel_coor = (-chirps/2:chirps/2)*radar_cube_struct.rfParams.dopplerResolutionMps;
angle_reso=180;
angle_coor = -angle_reso/2:angle_reso/2-1;
angle_coor_fft = asin((-90:1:90-1)/90)/pi*180;
%% range fft
figure(1);
% 相干叠加
range = squeeze(sum(radar_cube_all(1,:,:),2));
% range = squeeze(sum(sum(radar_cube_all,1),2));
% 非相干叠加
% range = squeeze(sum(sum(abs(radar_cube_all),1),2));
range_fft = 20*log2(squeeze(abs(range)));
% pick range peak
[range_max_val,range_max_idx] = max(range_fft);
%% doppler fft
subplot(2,1,2);
range_doppler_radarcube = fftshift(fft(radar_cube_all,[],1),1);
range_doppler_img = 20*log10(abs(squeeze(range_doppler_radarcube)));
range_doppler_img = squeeze(sum(range_doppler_img,2));
imagesc(range_coor,vel_coor,range_doppler_img);hold on
set(gca,'YDir','normal');
% pick velocity peak
vel_channel_data = abs(range_doppler_img(:,range_max_idx));
[vel_max_val, vel_max_idx] = max(vel_channel_data);
% do range doppler in each chirp*range curve
subplot(2,1,1);
for i=1:radar_cube_struct.rfParams.numDopplerBins
    range_curve = range_doppler_img(i,:);
    plot(range_coor, range_curve);hold on
    pause(0.01);    
end
