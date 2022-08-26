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
pack
c = 3e8; % m/s
setup_filename = './data/0805_angle/setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
%% 可调参数
range_cfar_thresh = 0.1;
doppler_cfar_thresh = 0.1;
frame_number = 1;
custome_peak_indx = 2;  % 调整这个值选取第几个峰为目标峰
MIMO = 1;
SKIP_read = 0;
%% 开始
frame = radar_cube_struct.data{frame_number};
if(MIMO)
    if(SKIP_read)
        radar_cube_tx1 = frame(1:3:end,:,:);
        radar_cube_tx2 = frame(2:3:end,:,:);
        radar_cube_tx3 = frame(3:3:end,:,:);
        radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);
    else
        radar_cube_tx1 = frame(1:128,:,:);
        radar_cube_tx2 = frame(129:256,:,:);
        radar_cube_tx3 = frame(257:end,:,:);
        radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);
    end
else
    radar_cube_all = frame;
end

[chirps,atennas,range_all_indx] = size(radar_cube_all);
range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(0:range_all_indx-1);
vel_coor = (-chirps/2:chirps/2-1)*radar_cube_struct.rfParams.dopplerResolutionMps;
angle_reso=180;
angle_coor = -angle_reso/2:angle_reso/2-1;
angle_coor_fft = asin((-90:1:90-1)/90)/pi*180;
%% 天线通道相位校准
load attena_calib_vec.mat
load attena_calib_ref.mat
attena_calib_vec_norm = attena_calib_ref./attena_calib_vec;
attena_calib_mat = repmat(attena_calib_vec_norm,[radar_cube_struct.rfParams.numDopplerBins,1,radar_cube_struct.rfParams.numRangeBins]);
radar_cube_all = radar_cube_all.*attena_calib_mat;

%% range fft
figure(1);
subplot(2,1,1);
% do range doppler in each chirp*range curve
range_doppler_img = squeeze(sum(abs(radar_cube_all),2));
bi_rangefft = zeros(1,256);
% 如果不把所有的chirp都sum起来，就这样操作
% for i=1:radar_cube_struct.rfParams.numDopplerBins
%     range_fft = 20*log10(range_doppler_img(i,:));
% 如果把所有的chirp都sum起来，就这样操作
for i=1:1
    range_fft = 20*log10(sum(range_doppler_img(i,:),1));
    bi_rangefft_temp = cfar_ca1D_square(range_fft,15,6,range_cfar_thresh,0);
    bi_rangefft = bi_rangefft+bi_rangefft_temp;
    plot(range_coor, range_fft);hold on

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
    
    pause(0.00001);
end

%% doppler fft
subplot(2,1,2);
range_doppler_radarcube = fftshift(fft(radar_cube_all,[],1),1);
range_doppler_img = 20*log10(abs(squeeze(range_doppler_radarcube)));
range_doppler_img = squeeze(sum(range_doppler_img,2));
imagesc(range_coor,vel_coor,range_doppler_img);hold on
set(gca,'YDir','normal');
% pick velocity peak
vel_channel_data = abs(range_doppler_img(:,peak_index_list));
[vel_max_val, vel_max_idx] = max(vel_channel_data);
plot(range_coor(peak_index_list),vel_coor(vel_max_idx), 'bo');
%% aoa
figure(3);
peak_number = length(peak_index_list);
row_num = 3;
col_num = ceil(peak_number/3);

for i = 1:peak_number
    TDM_Channel_data = zeros(1,angle_reso);
    range_peak_index = peak_index_list(i);
    vel_max_idx_single = vel_max_idx(i);
    TDM_Channel_data(1:12) = 20*log10(range_doppler_radarcube(vel_max_idx_single,:,range_peak_index));
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
figure(4);
subplot(1,2,1);
[Angle_coor_fft, Range_coor] = meshgrid(angle_coor_fft, range_coor);
mesh(Angle_coor_fft, Range_coor,AOA_heatmap_fft);view([0,0,1]);
% imagesc(angle_coor_fft, range_coor, AOA_heatmap_fft(1:end,:));
set(gca,'YDir','normal');
hold on;
subplot(1,2,2);
imagesc(angle_coor, range_coor, squeeze(AOA_heatmap));
set(gca,'YDir','normal');
hold on;

%% 静态测角
custome_peak_indx = peak_index_list(custome_peak_indx);
range_val = range_coor(custome_peak_indx);
% % fft
% angle_spectrum_fft = AOA_heatmap_fft(custome_peak_indx,:);
% [max_val,max_angle_index_fft] = max(angle_spectrum_fft);
% max_angle_fft = angle_coor_fft(max_angle_index_fft);
% subplot(1,2,1);plot3(max_angle_fft,range_val,max_val,'ro');
% % music
% angle_spectrum = AOA_heatmap(custome_peak_indx,:);
% [max_val,max_angle_index] = max(angle_spectrum);
% max_angle = angle_coor(max_angle_index);
% subplot(1,2,2);plot(max_angle,range_val,'ro');

%% 多通道相位校正


%% doppler compensation
% vel_bin = radar_cube_struct.rfParams.numDopplerBins;
% phi = 2*pi*(-vel_bin/2:vel_bin/2-1)/128/3;
% doppler_compensation_array = exp(-1i*pi*(-vel_bin/2:vel_bin/2-1)/128)';%多普勒相偏补偿
% doppler_compensation__mat = repmat(doppler_compensation_array,[1,12,radar_cube_struct.rfParams.numRangeBins]);
% radar_cube_all = radar_cube_all.*doppler_compensation__mat;
    

chirp_channel_mat = squeeze(radar_cube_all(:,:,custome_peak_indx));
% channel_mod = abs(chirp_channel_mat);
% channel_phase = angle(chirp_channel_mat)*180/pi;
% [chirps_number, channels_number] = size(chirp_channel_mat);
% for i = 1:chirps_number
%     subplot(1,2,1);
%     plot(channel_mod(i,:));title('幅值分布');hold on
%     subplot(1,2,2);
%     plot(channel_phase(i,:));title('相位分布');hold on
%     pause(0.01);
% end
% 
% attena_calib_vec = mean(chirp_channel_mat,1);
% attena_calib_ref = mean(chirp_channel_mat(:,1),1);
% 
% attena_calib_vec_norm = attena_calib_ref./attena_calib_vec;
% attena_calib_mat = repmat(attena_calib_vec_norm,[chirps_number,1,radar_cube_struct.rfParams.numRangeBins]);
% % save attena_calib_vec.mat attena_calib_vec
% % save attena_calib_ref.mat attena_calib_ref
% 
% radar_cube_all_calib = radar_cube_all.*attena_calib_mat;
% 
% 
% chirp_channel_mat = squeeze(radar_cube_all_calib(:,:,custome_peak_indx));
% channel_mod = abs(chirp_channel_mat);
% channel_phase = angle(chirp_channel_mat)*180/pi;
% [chirps_number, channels_number] = size(chirp_channel_mat);
% for i = 20:chirps_number
%     subplot(1,2,1);
%     plot(channel_mod(i,:));title('幅值分布');hold on
%     subplot(1,2,2);
%     plot(channel_phase(i,:));title('相位分布');hold on
%     pause(0.01);
% end
% 
% hold off;

figure(5);
for i=1:128
    signal1 = squeeze(chirp_channel_mat(i,1:4));
    signal2 = squeeze(chirp_channel_mat(i,5:8));
    signal3 = squeeze(chirp_channel_mat(i,9:12));
    
%     signal_random = squeeze(chirp_channel_mat(i,4:));
    
    signal_set1 = [signal1,signal3];
    signal_set2 = [signal1,signal2];
    signal_set3 = [signal2,signal3];
    signal_set4 = [signal2,signal1];
    signal_set5 = [signal3,signal2];
    signal_set6 = [signal3,signal1];
    signal_set7 = [signal1,signal2,signal3];
    subplot(4,2,1);plot(musicAlg(signal1',angle_reso));hold on;title('单组')
    subplot(4,2,2);plot(musicAlg(signal_set1',angle_reso));hold on;title('1+3组')
    subplot(4,2,3);plot(musicAlg(signal_set2',angle_reso));hold on;title('1+2组')
    subplot(4,2,4);plot(musicAlg(signal_set3',angle_reso));hold on;title('2+3组')
    subplot(4,2,5);plot(musicAlg(signal_set4',angle_reso));hold on;title('2+1组')
    subplot(4,2,6);plot(musicAlg(signal_set5',angle_reso));hold on;title('3+2组')
    subplot(4,2,7);plot(musicAlg(signal_set6',angle_reso));hold on;title('3+1组')
    subplot(4,2,8);plot(musicAlg(signal_set7',angle_reso));hold on;title('全组')
    pause(0.0001);
    
end
hold off
