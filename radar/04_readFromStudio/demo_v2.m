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
% update: 2022-08-09，demo_v2是用来读取动态数据，并可视化的
%
%
close all
clear all
c = 3e8; % m/s
setup_filename = './data/0809/setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
%% 可调参数
range_cfar_thresh = 0.1;
doppler_cfar_thresh = 0.1;
frame_number = 3;
custome_peak_indx = 1;  % 调整这个值选取第几个峰为目标峰
%% 开始
figure(1);
for frame_number=1:radar_cube_struct.dim.numFrames
    frame = radar_cube_struct.data{frame_number};
    radar_cube_tx1 = frame(1:3:end,:,:);
    radar_cube_tx2 = frame(2:3:end,:,:);
    radar_cube_tx3 = frame(3:3:end,:,:);
    radar_cube_all = cat(2,radar_cube_tx1, radar_cube_tx2, radar_cube_tx3);

    [chirps,atennas,range_all_indx] = size(radar_cube_all);
    range_coor = c/(2*radar_cube_struct.rfParams.bandwidth*1e9)*(0:range_all_indx-1);
    vel_coor = (-chirps/2:chirps/2)*radar_cube_struct.rfParams.dopplerResolutionMps;
    angle_reso=180;
    angle_coor = -angle_reso/2:angle_reso/2-1;
%     %% 天线通道相位校准
%     load attena_calib_vec.mat
%     load attena_calib_ref.mat
%     attena_calib_vec_norm = attena_calib_ref./attena_calib_vec;
%     attena_calib_mat = repmat(attena_calib_vec_norm,[radar_cube_struct.rfParams.numDopplerBins,1,radar_cube_struct.rfParams.numRangeBins]);
%     radar_cube_all = radar_cube_all.*attena_calib_mat;

    %% range fft
    subplot(2,3,[1,2]);
    range = squeeze(sum(radar_cube_all(1,:,:),2));
    range_fft = 20*log10(squeeze(abs(range)));
    % pick range peak
    [range_max_val,range_max_idx] = max(range_fft);
    plot(range_coor,range_fft, 'b');hold on
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

    plot(range_coor(peak_index_list),peak_val_list, 'rp');hold on
    %% doppler fft
    subplot(2,3,3);
    range_doppler_radarcube = fftshift(fft(radar_cube_all,[],1),1);
    range_doppler_img = 20*log10(abs(squeeze(range_doppler_radarcube)));
    range_doppler_img = squeeze(sum(range_doppler_img,2));
    imagesc(range_coor,vel_coor,range_doppler_img);hold on
    set(gca,'YDir','normal');
    
    %% aoa
    AOA_heatmap = zeros(radar_cube_struct.rfParams.numRangeBins, angle_reso);
    range_doppler_peak_mask = zeros(radar_cube_struct.rfParams.numRangeBins, radar_cube_struct.rfParams.numDopplerBins);
    for i = 1:radar_cube_struct.rfParams.numRangeBins
        bi_ = bi_rangefft(i);
        if bi_ > 0
            doppler_peak = range_doppler_img(:,i);
            doppler_peak = cfar_ca1D_square(doppler_peak,12,6,doppler_cfar_thresh,0);
            bi_doppler_peak = doppler_peak == 1;
            range_doppler_peak_mask(i,:) = doppler_peak;
            aoa_signal = radar_cube_all(bi_doppler_peak,:,i);
            aoa_result = musicAlg(aoa_signal',angle_reso);
            AOA_heatmap(i,:) = aoa_result;
        end
    end
    subplot(2,3,6);imagesc(range_coor,vel_coor,range_doppler_peak_mask');set(gca,'YDir','normal');
    new_radar_cube_all = zeros(radar_cube_struct.rfParams.numDopplerBins,angle_reso,radar_cube_struct.rfParams.numRangeBins);
    new_radar_cube_all(:,1:12,:) = radar_cube_all;
    AOA_heatmap_fft = squeeze(20*log10(sum(abs(fft(new_radar_cube_all,[],2)),1)));
    AOA_heatmap_fft = fliplr(AOA_heatmap_fft');

    subplot(2,3,4);
    imagesc(angle_coor, range_coor, AOA_heatmap_fft);
    set(gca,'YDir','normal');
    subplot(2,3,5);
    imagesc(angle_coor, range_coor, squeeze(AOA_heatmap));
    set(gca,'YDir','normal');
    %% 输出最终检测点结果
    range_list = [];
    angle_list_fft = [];
    angle_list_music = [];
    x_fft = [];
    y_fft = [];
    x_music = [];
    y_music = [];
    for i = 1:length(peak_index_list)
        range_index = peak_index_list(i);
        range=range_coor(range_index);
        range_list = [range_list, range];
        %% fft aoa result
        [max_aoa_val_fft,max_aoa_index_fft] = max(AOA_heatmap_fft(range_index,:));
        angle_fft = angle_coor(max_aoa_index_fft);
        angle_list_fft = [angle_list_fft, angle_fft];
        x_fft = [x_fft, range*sin(angle_fft*pi/180)];
        y_fft = [y_fft, range*cos(angle_fft*pi/180)];
        %% music aoa result
        [max_aoa_val_music,max_aoa_index_music] = max(AOA_heatmap(range_index,:));
        angle_music = angle_coor(max_aoa_index_music);
        angle_list_music = [angle_list_music, angle_music];
        x_music = [x_music, range*sin(angle_music*pi/180)];
        y_music = [y_music, range*cos(angle_music*pi/180)];
    end
    figure(2);
    subplot(1,2,1);
    scatter(x_fft, y_fft,'b','filled');hold on;
    scatter(x_music, y_music,'g','filled');hold on;
    legend('fft','music');
    subplot(1,2,2);
    scatter(angle_list_fft, range_list,'b','filled');hold on;
    scatter(angle_list_music, range_list,'g','filled');hold on;
    legend('fft','music');
    pause(0.01);
    clf(figure(1));

    %% 静态测角
%     custome_peak_indx = peak_index_list(custome_peak_indx);
%     range_val = range_coor(custome_peak_indx);
%     % fft
%     angle_spectrum_fft = AOA_heatmap_fft(custome_peak_indx,:);
%     [max_val,max_angle_index_fft] = max(angle_spectrum_fft);
%     max_angle_fft = angle_coor(max_angle_index_fft);
%     figure(2);subplot(1,2,1);plot(max_angle_fft,range_val,'ro');hold on
%     % music
%     angle_spectrum = AOA_heatmap(custome_peak_indx,:);
%     [max_val,max_angle_index] = max(angle_spectrum);
%     max_angle = angle_coor(max_angle_index);
%     figure(2);subplot(1,2,1);plot(max_angle,range_val,'ro');hold on

end

