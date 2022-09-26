%
% 简介：
% 利用ti给的rawDataReader.m
% 来进行studio数据的读取
% 
% 本脚本的功能是进行动态测角谱分析
% 逻辑是，先读取数据，进行数据组包
% 然后进行每一个chirp下的距离FFT，doppler fft和角度music分析
% 与静态角度分析不同的是，并不是每一个chirp都是我们需要的测角chirp了，我们需要取出特定的chirp进行测角music分析
% 这样的工作流程就已经是最接近真实情况下的毫米波雷达工作流程了，因此可能存在的问题有：测角的doppler补偿问题
% 2022-09-26

close all
clear all
c = 3e8; % m/s
setup_filename = './data/0923_高级仪器测试/studio配置/1.setup.json';
radar_cube_struct = rawDataReader(setup_filename, 'temp_raw_bin.bin', 'temp_cube', 0);
%% 可调参数
range_cfar_thresh = 0.1;
doppler_cfar_thresh = 0.1;
angle_cfar_thresh = 0.5;
frame_number = 1;
custome_peak_indx = 1;  % 调整这个值选取第几个峰为目标峰
MIMO = 1;
SKIP_read = 1;
MIMO_2attena = 1;   % 仅仅使用TX0&TX2
ANGLE_RESO = 180;
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
figure(1);subplot(2,2,1);
imagesc(range_coor,vel_coor,log2(range_doppler_img));set(gca,'YDir','normal');title('range doppler heatmap');
% 如果不把所有的chirp都sum起来，就这样操作
% for i=1:radar_cube_struct.rfParams.numDopplerBins
%     range_fft = 20*log10(range_doppler_img(i,:));
% 如果把所有的chirp都sum起来，就这样操作
%% ----------------------------------------
% 在range doppler图上进行2D CFAR
% 选择进行2次1D CFAR的理由是，这样的cfar比2D CFAR更具有鲁棒性
% ----------------------------------------
bi_range_doppler_img_range = zeros(size(range_doppler_img));
bi_range_doppler_img_doppler = zeros(size(range_doppler_img));
figure(1);
% range dimension cfar
for i=1:chirps
    range_fft = 20*log10(range_doppler_img(i,:));
    bi_rangefft_temp = cfar_ca1D_square(range_fft,15,6,range_cfar_thresh,0);
    bi_range_doppler_img_range(i,:) = bi_rangefft_temp;
end
% doppler dimension cfar
for i=1:range_all_indx
    doppler_fft = 20*log10(range_doppler_img(:,i));
    bi_rangefft_temp = cfar_ca1D_square(doppler_fft,15,6,doppler_cfar_thresh,0);
    bi_range_doppler_img_doppler(:,i) = bi_rangefft_temp;
end
bi_range_doppler_img = bi_range_doppler_img_range&bi_range_doppler_img_doppler;
subplot(2,2,3);imagesc(range_coor,vel_coor,bi_range_doppler_img);set(gca,'YDir','normal');hold on;title('range doppler peak map');
[doppler_row,range_col] = find(bi_range_doppler_img);
% subplot(2,2,3);plot(range_coor(range_col),vel_coor(doppler_row),'rp');hold on;plot(range_coor(range_col),vel_coor(doppler_row));
%% ----------------------------------------
% 对在range-doppler图上取到的动态目标进行角度分析
% ----------------------------------------
angle_coor = -ANGLE_RESO/2:ANGLE_RESO/2-1;
% 遍历每个距离，在特定距离下进行music测角，构造radar range-azimuth heatmap
range_azimuth_heatmap = zeros(range_all_indx, ANGLE_RESO);
bi_range_azimuth_heatmap = zeros(range_all_indx, ANGLE_RESO);
angle_spectrum = zeros(1, ANGLE_RESO);
num_index = length(doppler_row);
old_bi_range_doppler_img_val = 0;
music_input_mat = [];
% result list
range_list = [];
max_angle_spectrum_peak_value_list = [];
max_angle_spectrum_peak_angle_value_list = [];
for i=1:range_all_indx
    for j=1:chirps
        bi_range_doppler_img_val = bi_range_doppler_img(j,i);
        delta_doppler_index = bi_range_doppler_img_val-old_bi_range_doppler_img_val;
        old_bi_range_doppler_img_val = bi_range_doppler_img_val;
        if delta_doppler_index == 1 && bi_range_doppler_img_val == 1
            channel_signal = radar_cube_all(j,:,i);
            music_input_mat = [music_input_mat; channel_signal];
        elseif delta_doppler_index == -1 && bi_range_doppler_img_val == 0
            % do music here
%             angle_spectrum = angle_spectrum+musicAlg(music_input_mat', ANGLE_RESO);
            range_azimuth_heatmap(i,:) = range_azimuth_heatmap(i,:)+musicAlg(music_input_mat', ANGLE_RESO);
            music_input_mat = [];
        end
    end
    % ---------------------------------------
    % CFAR在range-azimuth heatmap上进行目标提取
    % ---------------------------------------
    bi_range_azimuth_heatmap(i,:) = cfar_ca1D_square(range_azimuth_heatmap(i,:),15,6,angle_cfar_thresh,0);
    % peak prune
    max_angle_spectrum_peak_value = 0;
    old_bi_angle_spectrum_value = 0;
    for k=1:ANGLE_RESO
        bi_angle_spectrum_value = bi_range_azimuth_heatmap(i,k);
        delta_bi_angle_spectrum_value = bi_angle_spectrum_value-old_bi_angle_spectrum_value;
        old_bi_angle_spectrum_value = bi_angle_spectrum_value;
        if bi_angle_spectrum_value==1
%         if delta_bi_angle_spectrum_value==1 && bi_angle_spectrum_value==1
            angle_spectrum_value = range_azimuth_heatmap(i,k);
            if angle_spectrum_value > max_angle_spectrum_peak_value
                max_angle_spectrum_peak_value = angle_spectrum_value;
                max_angle_spectrum_peak_angle_value = angle_coor(k);
            end
        elseif delta_bi_angle_spectrum_value==-1 && bi_angle_spectrum_value==0
            % 完成了一个峰所有的值的检索，这时候剩下的max peak就是peak prune的结果
            range_list = [range_list, range_coor(i)];
            max_angle_spectrum_peak_value_list = [max_angle_spectrum_peak_value_list, max_angle_spectrum_peak_value];
            max_angle_spectrum_peak_angle_value_list = [max_angle_spectrum_peak_angle_value_list, max_angle_spectrum_peak_angle_value];
            max_angle_spectrum_peak_value = 0;
        end
    end
end
subplot(2,2,[2,4]);
imagesc(angle_coor,range_coor,range_azimuth_heatmap);set(gca,'YDir','normal');hold on;title('range azimuth heatmap');
[range_row,angle_col] = find(bi_range_azimuth_heatmap);
% subplot(2,2,[2,4]);
% plot(angle_coor(angle_col),range_coor(range_row),'gp');hold on
subplot(2,2,[2,4]);
plot(max_angle_spectrum_peak_angle_value_list,range_list,'rp');hold on