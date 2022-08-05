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
setup_filename = './data/0805_angle/setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
%% 可调参数
range_cfar_thresh = 0.25;
doppler_cfar_thresh = 0.1;
frame_number = 1;
custome_peak_indx = 2;  % 调整这个值选取第几个峰为目标峰
%% 开始
frame = radar_cube_struct.data{frame_number};
radar_cube_tx1 = frame(1:3:end,:,:);
radar_cube_tx2 = frame(2:3:end,:,:);
radar_cube_tx3 = frame(3:3:end,:,:);
radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);
figure(1);
[chirps,atennas,range_all_indx] = size(radar_cube_all);
range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(0:range_all_indx-1);
vel_coor = (-chirps/2:chirps/2)*radar_cube_struct.rfParams.dopplerResolutionMps;
angle_reso=180;
angle_coor = -angle_reso/2:angle_reso/2-1;
%% range fft
subplot(2,1,1);
range = squeeze(sum(radar_cube_all(1,:,:),2));
range_fft = 20*log10(squeeze(abs(range)));
% pick range peak
[range_max_val,range_max_idx] = max(range_fft);
plot(range_coor,range_fft, 'r');hold on
% use cfar to pick peak
[bi_rangefft] = cfar_ca1D_square(range_fft,15,6,doppler_cfar_thresh,0);
% use peak prune to pick a single peak index
peak_index_list = [];
peak_val_list = [];
peak_val = 0;
peak_index = 0;
for i = 1:length(bi_rangefft)
    bi_ = bi_rangefft(i);
    if bi_ > 0
        val_rangefft = range_fft(i);
        if peak_val < val_rangefft
            peak_index = i;
            peak_val = val_rangefft;
        end
    else
        if peak_index == 0
            continue;
        end
        peak_index_list = [peak_index_list, peak_index];
        peak_val_list = [peak_val_list, peak_val];
        peak_val = 0;
        peak_index = 0;
    end
end

plot(range_coor(peak_index_list),peak_val_list, 'bo');hold on
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
plot(range_coor(peak_index_list),vel_coor(vel_max_idx), 'bo');
%% aoa
figure(2);
peak_number = length(peak_index_list);
row_num = 3;
col_num = ceil(peak_number/3);

for i = 1:peak_number
    TDM_Channel_data = zeros(1,angle_reso);
    range_peak_index = peak_index_list(i);
    TDM_Channel_data(1:12) = 20*log10(range_doppler_radarcube(vel_max_idx,:,range_peak_index));
%     aoa_spectrum_analysis = musicAlg(TDM_Channel_data',angle_reso);
    aoa_spectrum_analysis = abs(fftshift(fft(TDM_Channel_data)));
    subplot(row_num,col_num*2,2*i-1);
    plot(abs(TDM_Channel_data));
    subplot(row_num,col_num*2,2*i);
%     plot(anlge_coor, aoa_spectrum_analysis);
    plot([-angle_reso/2:angle_reso/2-1]*180/36,abs(aoa_spectrum_analysis));

end

AOA_heatmap = zeros(radar_cube_struct.rfParams.numRangeBins, angle_reso);
for i = 1:radar_cube_struct.rfParams.numRangeBins
    bi_ = bi_rangefft(i);
    if bi_ > 0
        doppler_peak = range_doppler_img(:,i);
        bi_doppler_peak = cfar_ca1D_square(doppler_peak,12,6,doppler_cfar_thresh,0)==1;
        aoa_signal = radar_cube_all(bi_doppler_peak,:,i);
        aoa_result = musicAlg(aoa_signal',angle_reso);
        AOA_heatmap(i,:) = aoa_result;        
    end
end

new_radar_cube_all = zeros(radar_cube_struct.rfParams.numDopplerBins,angle_reso,radar_cube_struct.rfParams.numRangeBins);
new_radar_cube_all(:,1:12,:) = radar_cube_all;
AOA_heatmap_fft = squeeze(20*log10(sum(abs(fft(new_radar_cube_all,[],2)),1)));
AOA_heatmap_fft = fliplr(AOA_heatmap_fft');
figure(3);
subplot(1,2,1);
imagesc(angle_coor, range_coor, AOA_heatmap_fft);
set(gca,'YDir','normal');
hold on;
subplot(1,2,2);
imagesc(angle_coor, range_coor, squeeze(AOA_heatmap));
set(gca,'YDir','normal');
hold on;

%% 静态测角
custome_peak_indx = peak_index_list(custome_peak_indx);
range_val = range_coor(custome_peak_indx);
% fft
angle_spectrum_fft = AOA_heatmap_fft(custome_peak_indx,:);
[max_val,max_angle_index_fft] = max(angle_spectrum_fft);
max_angle_fft = angle_coor(max_angle_index_fft);
subplot(1,2,1);plot(max_angle_fft,range_val,'ro');
% music
angle_spectrum = AOA_heatmap(custome_peak_indx,:);
[max_val,max_angle_index] = max(angle_spectrum);
max_angle = angle_coor(max_angle_index);
subplot(1,2,2);plot(max_angle,range_val,'ro');


%% 静态误差曲线
gt_angle = [-45,-30:5:-5,5:5:30];
fft_angle = [-44,-46,-46,-47,-49,-50,89,89,87,87,46,45,44];
music_angle = [-32,-30,-30,-30,-2,-1,1,1,1,1,28,29,29];
figure(4);
subplot(2,1,1);
plot(gt_angle,gt_angle,'ro-');hold on
plot(gt_angle,fft_angle,'go-');hold on
plot(gt_angle,music_angle,'bo-');hold on
legend('gt','fft','music');
subplot(2,1,2);
% plot(fft_angle-gt_angle,'go-');hold on
plot(gt_angle,music_angle-gt_angle,'bo-');hold on
legend('music');
hold off;
