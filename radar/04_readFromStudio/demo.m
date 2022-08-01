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
% 2. 保存raw_adc.bin与xx.setup.json文件。其中前者为数据文件，后者为当前数据的配置文件
% 2.1 注意修改setup.json中的configUsed和fileBasePath到数据对应目录下，数据名称processedFileName修改保存名称
% 3. 将数据与配置文件放入data文件夹，进行数据可视化
%
close all
clear all
c = 3e8; % m/s
setup_filename = './data/1.setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
frame = radar_cube_struct.data{1};
radar_cube_tx1 = frame(1:3:96,:,:);
radar_cube_tx2 = frame(2:3:96,:,:);
radar_cube_tx3 = frame(3:3:96,:,:);
radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);
figure(1);
[chirps,atennas,range_all_indx] = size(radar_cube_all);
range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(1:range_all_indx);
vel_coor = (-chirps/2:chirps/2)*radar_cube_struct.rfParams.dopplerResolutionMps;
angle_reso=180;
anlge_coor = -angle_reso/2:angle_reso/2-1;
%% range fft
subplot(3,1,1);
range = radar_cube_all(1,1,:);
range_fft = 20*log10(squeeze(abs(range)));
plot(range_coor,range_fft, 'r');hold on
% pick range peak
[range_max_val,range_max_idx] = max(range_fft);
obj_range = range_coor(range_max_idx);
plot(obj_range,range_max_val, 'bo');
% legend(obj_range,max_val);
%% doppler fft
subplot(3,1,2);
range_doppler_radarcube = fftshift(fft(radar_cube_all,[],1),1);
range_doppler_img = 20*log10(abs(squeeze(range_doppler_radarcube(:,1,:))));
imagesc(range_coor,vel_coor,range_doppler_img);hold on
set(gca,'YDir','normal');
% pick velocity peak
vel_channel_data = abs(range_doppler_img(:,range_max_idx));
[vel_max_val, vel_max_idx] = max(vel_channel_data);
plot(obj_range,vel_coor(vel_max_idx), 'bo');
%% aoa
subplot(3,1,3);
TDM_Channel_data = 20*log10(radar_cube_all(vel_max_idx,:,range_max_idx));
aoa_spectrum_analysis = musicAlg(TDM_Channel_data',angle_reso);
plot(anlge_coor, aoa_spectrum_analysis);

hold off;
