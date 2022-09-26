%
% 简介：
% 利用ti给的rawDataReader.m
% 来进行studio数据的读取
% 
% 本脚本的功能是进行静态测角谱分析
% 逻辑是，先读取数据，进行数据组包
% 然后进行每一个chirp下的距离FFT，doppler fft和角度music分析
% 
% 
% 
close all
clear all
c = 3e8; % m/s
setup_filename = './data/0923_高级仪器测试/studio配置/1.setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
%% 可调参数
range_cfar_thresh = 0.02;
doppler_cfar_thresh = 0.1;
frame_number = 1;
custome_peak_indx = 1;  % 调整这个值选取第几个峰为目标峰
MIMO = 1;
SKIP_read = 1;
MIMO_2attena = 1;   % 仅仅使用TX0&TX2
%% 组包数据
frame = radar_cube_struct.data{frame_number};
if(MIMO)
    if(MIMO_2attena)
        if(SKIP_read)
            radar_cube_tx1 = frame(1:2:end,:,:);
            radar_cube_tx2 = frame(2:2:end,:,:);
            radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2);
        else
            radar_cube_tx1 = frame(1:128,:,:);
            radar_cube_tx2 = frame(129:256,:,:);
            radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2);
        end
    else
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
    end
else
    radar_cube_all = frame;
end
%% range fft
[chirps,atennas,range_all_indx] = size(radar_cube_all);
range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(0:range_all_indx-1);


% do range doppler in each chirp*range curve
% do doppler fft
radar_cube_all = fftshift(fft(radar_cube_all,[],1),1);
% radar_cube_all = fft(radar_cube_all,[],1);
range_doppler_img = squeeze(sum(abs(radar_cube_all),2));
vel_coor = (-chirps/2:chirps/2)*radar_cube_struct.rfParams.dopplerResolutionMps;
figure(10);imagesc(range_coor,vel_coor,log2(range_doppler_img));set(gca,'YDir','normal');
bi_rangefft = zeros(1,256);
% 如果不把所有的chirp都sum起来，就这样操作
% for i=1:radar_cube_struct.rfParams.numDopplerBins
%     range_fft = 20*log10(range_doppler_img(i,:));
% 如果把所有的chirp都sum起来，就这样操作
figure(1);
for i=1:128
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
    if MIMO
        subplot(5,2,[1,2]);
    else
        subplot(3,2,[1,2]);
    end
    plot(range_coor(peak_index_list),peak_val_list, 'bo');hold on
    
    pause(0.00001);
end
%% 天线通道相位校准
% load attena_calib_vec.mat
% load attena_calib_ref.mat
% attena_calib_vec_norm = attena_calib_ref./attena_calib_vec;
% attena_calib_mat = repmat(attena_calib_vec_norm,[radar_cube_struct.rfParams.numDopplerBins,1,radar_cube_struct.rfParams.numRangeBins]);
% radar_cube_all = radar_cube_all.*attena_calib_mat;


%% 计算幅值相位custome_peak_indx = peak_index_list(custome_peak_indx);
range_val = range_coor(custome_peak_indx);
% custome_peak_indx = 8;
chirp_channel_mat = squeeze(radar_cube_all(:,:,custome_peak_indx));
channel_mod = abs(chirp_channel_mat);
channel_phase = angle(chirp_channel_mat)*180/pi;
[chirps_number, channels_number] = size(chirp_channel_mat);
for i = 1:chirps_number
    if MIMO
        subplot(5,2,3);
    else
        subplot(3,2,3);
    end
    plot(channel_mod(i,:));title('幅值分布');hold on
    if MIMO
        subplot(5,2,4);
    else
        subplot(3,2,4);
    end
    plot(channel_phase(i,:));title('相位分布');hold on
    pause(0.001);
end
%% 计算天线补偿值
% % 求天线复信号均值
% attena_calib_vec = mean(chirp_channel_mat,1);
% % 选择最后一根天线为参考
% attena_calib_ref = mean(chirp_channel_mat(:,12),1);
% % 计算相位补偿信号：参考/每根
% attena_calib_vec_norm = attena_calib_ref./attena_calib_vec;
% attena_calib_mat = repmat(attena_calib_vec_norm,[chirps_number,1,radar_cube_struct.rfParams.numRangeBins]);
% % 保存天线相位补偿矩阵
% % save attena_calib_vec.mat attena_calib_vec
% % save attena_calib_ref.mat attena_calib_ref
% % 使用天线相位补偿矩阵
% radar_cube_all_calib = radar_cube_all.*attena_calib_mat;
%% 再画校准后的结果
% chirp_channel_mat = squeeze(radar_cube_all_calib(:,:,custome_peak_indx));
% channel_mod = abs(chirp_channel_mat);
% channel_phase = angle(chirp_channel_mat)*180/pi;
% [chirps_number, channels_number] = size(chirp_channel_mat);
% for i = 1:chirps_number
%     subplot(1,2,1);
%     plot(channel_mod(i,:));title('幅值分布');hold on
%     subplot(1,2,2);
%     plot(channel_phase(i,:));title('相位分布');hold on
%     pause(0.001);
% end

hold off;

%% 计算角度
angle_reso=180;
for i=1:128
    if MIMO
        if MIMO_2attena
            signal1 = squeeze(chirp_channel_mat(i,1:4));
            signal2 = squeeze(chirp_channel_mat(i,5:8));
            signal_set2 = [signal1,signal2];
%             subplot(5,2,[5,6,7,8]);plot(musicAlg(signal1',angle_reso));hold on;title('单组')
            subplot(5,2,[5,6,7,8]);plot(musicAlg(signal_set2',angle_reso));hold on;title('0+2组')
%             subplot(5,2,[9,10]);plot(musicAlg(signal_set7',angle_reso));hold on;title('全组')
        else
            signal1 = squeeze(chirp_channel_mat(i,1:4));
            signal2 = squeeze(chirp_channel_mat(i,5:8));
            signal3 = squeeze(chirp_channel_mat(i,9:12));
            signal_set1 = [signal1,signal3];
            signal_set2 = [signal1,signal2];
            signal_set3 = [signal2,signal3];
            signal_set7 = [signal1,signal2,signal3];
            subplot(5,2,[5,6]);plot(musicAlg(signal_set1',angle_reso));hold on;title('1+3组')
            subplot(5,2,[7,8]);plot(musicAlg(signal_set2',angle_reso));hold on;title('1+2组')
            subplot(5,2,9);plot(musicAlg(signal2',angle_reso));hold on;title('单组')
            subplot(5,2,10);plot(musicAlg(signal_set7',angle_reso));hold on;title('全组')
        end
        pause(0.0001);
    else 
        signal1 = squeeze(chirp_channel_mat(i,:));
        subplot(3,2,[5,6]);plot(musicAlg(signal1',angle_reso));hold on;title('单组');
    end
    
    
end
hold off