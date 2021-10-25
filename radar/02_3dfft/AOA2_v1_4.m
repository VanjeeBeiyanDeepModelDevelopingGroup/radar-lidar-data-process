%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     AOA algorithm                                                                            %%%
%%%     dopplerOut- doppler FFT的结果，作为函数的输入                                            %%%  
%%%                 矩阵维度分别为 numADCSamples, dopplerBin_num, RX_num, TX_num, frame_num      %%% 
%%%     CFAROut- CFAR后得到峰值的每帧行列数                                                      %%%
%%%              一共帧数个元胞数组，每个元胞数组中的第一列是range索引，第二列为doppler索引      %%%
%%%     distanceCoor- 距离坐标                                                                   %%%
%%%     TX_num- 发射天线数量                                                                     %%%
%%%     RX_num- 接收天线数量                                                                     %%%
%%%     numADCSamples- range bin数量                                                             %%%
%%%     dopplerBin_num- doppler FFT bin容量                                                      %%%
%%%     angleBin_num- 角度估计FFT bin容量                                                        %%%
%%%     MAX_VEL_ENH_PROCESSING- 是否使能MAX VEL功能的标志位                                      %%%
%%%     azimuthOut- 水平角度3D FFT结果                                                           %%%
%%%     elevOut- 垂直角度3D FFT结果                                                              %%%
%%%                                                                                              %%%
%%%     Created by 李嘉宝 2021.03.18 version 1.4                                                 %%%
%%%     修改部分： 去除cfar结果筛选，对整个数据进行angle fft                                     %%%
%%%                画range-angle heat map                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ azimuthOut, elevOut, pcStrc ] = AOA2_v1_4( dopplerOut, CFAROut, distanceCoor,velocityCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frameIndex,lidarDataFrame,lineId)
% function [ azimuthOut, elevOut, agl_grid, ] = AOA2_v1_4( dopplerOut, CFAROut, distanceCoor, TX_num, RX_num, numADCSamples, dopplerBin_num, angleBin_num, MAX_VEL_ENH_PROCESSING, frameIndex,lidarDataFrame)
    azimuthOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
    elevOut = zeros(numADCSamples, angleBin_num, dopplerBin_num);
%     azimuthOut = zeros(numADCSamples, 512, dopplerBin_num);
%     elevOut = zeros(numADCSamples, 512, dopplerBin_num);
%      figure(12);
     for r = 1: numADCSamples
         
         for c = 1: dopplerBin_num

             if TX_num == 3

                 Temp1 = zeros(1, angleBin_num);
                 Temp2 = zeros(1, angleBin_num);

                 Temp3 = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp3 = RxPhaseBiasCompensation(Temp3, TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp3 = dopplerCompensation(Temp3, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);
                 Temp1(1: ((TX_num - 1) * RX_num)) = Temp3(1: ((TX_num - 1) * RX_num)); % 8个值
                 Temp2(1: RX_num) = Temp3((((TX_num - 1) * RX_num) + 1): (TX_num * RX_num)); % 4个值
                 %% 加窗
%                  angleWin = hanning(angleBin_num);  % 32个chirp
%                  angleWin = angleWin(1: (angleBin_num / 2));
%                  angleWinLen               = length(angleWin);
%                  angleWindowCoeffVec       = ones(angleBin_num, 1);
%                  angleWindowCoeffVec(1:angleWinLen) = angleWin;
%                  angleWindowCoeffVec(angleBin_num-angleWinLen+1:angleBin_num) = angleWindowCoeffVec(angleWinLen:-1:1);
%                  angleWin = angleWindowCoeffVec;
%                  Temp1 = bsxfun(@times, Temp1, angleWin.');
                 azimuthOut(r, :, c) = fftshift(fft(Temp1));
%                  [music,f1] = pmusic(Temp1,32,angleBin_num,'whole');
%                  azimuthOut(r, :, c) = 20*log10(abs(music));
%                  plot(abs(azimuthOut(r, :, c)));hold on
%                  azimuthOut(r, :, c) = fft(Temp1);
                 elevOut(r, :, c) = fftshift(fft(Temp2));
%                  [elevOut(r, :, c),f2] = pmusic(Temp2,32,angleBin_num,'whole');

             else

                 Temp1 = zeros(1, angleBin_num);
                 Temp1(1: TX_num * RX_num) = reshape(dopplerOut(r, c, :, :), [1, TX_num * RX_num]);
                 Temp1(1: (TX_num * RX_num)) = RxPhaseBiasCompensation(Temp1(1: (TX_num * RX_num)), TX_num, RX_num, MAX_VEL_ENH_PROCESSING);
                 Temp1 = dopplerCompensation(Temp1, c, TX_num, RX_num, dopplerBin_num, MAX_VEL_ENH_PROCESSING);

                 azimuthOut(r, :, c) = fftshift(fft(Temp1));
%                 [azimuthOut(r, :, c),f1] = pmusic(Temp1,32,angleBin_num,'whole');

             end
         
         end

     end
     
     w = linspace(-1,1,angleBin_num); % angle_grid
     agl_grid = asin(w)*180/pi+5; % [-1,1]->[-pi/2,pi/2]
     noisefloor = db2pow(-80);
     Xpow = abs(azimuthOut(:, :, :)).^2;
     %% 画3d cube
%      Xpow_show = 20*log10(Xpow);
%      maxXpow = max(Xpow_show(:));
%      minXpow = min(Xpow_show(:));
%      XpowNorm = (Xpow_show-minXpow)/(maxXpow-minXpow);
%      figure(4);
%      [range,angle,doppler] = size(Xpow_show);
%      [r,a,d] = ndgrid(1:range,1:angle,1:doppler);
%      s = 10*(1-XpowNorm(:))+1;
%      scatter3(r(:),a(:),d(:),s(:),1-XpowNorm(:));
%      xlim([1,512]);
%      ylim([1,64]);
%      zlim([1,32]);
%      colormap pink
%      set(gca,'color','b');
    %%
%      Xpow = azimuthOut;
     % 画出每个速度下的512*64距离*角度图
     figure(3);
     for i=1:dopplerBin_num
         subplot(4,8,i);
%          fprintf("doppler: %d\n",i);
         imagesc(Xpow(:,:,i));
         titlename = num2str(velocityCoor(i));
         title(titlename);
         axis off
         set(gca,'YDir','normal');
%          pause(1);
     end
     
%      figure(9);
% %      newXpow = [];
% %      num = 32;
%      for i=1:3
%          subplot(1,3,i);
%          imagesc(Xpow(:,:,i+23));
% %          newXpow(:,:,i) = Xpow(:,:,i+23);
%          titlename = num2str(velocityCoor(i+23));
%          title(titlename);
%          axis off
%          set(gca,'YDir','normal');
%      end
%      newXpow_sum = squeeze(sum(newXpow,3)/size(newXpow,3));
%      Xpow_max = max(Xpow,[],3); % 这里跟我提供的代码不一样，他里面是sum起来，但是取max似乎是差不多的意义
%      Xpow_sum = squeeze(sum(Xpow,3)/size(Xpow,3));
%      Xpow_max_snr = pow2db(Xpow_max/noisefloor);
%      Xpow_sum_snr = pow2db(Xpow_sum/noisefloor);
%      newNoiseFloor = db2pow(-40);
%      newXpow_sum_snr = pow2db(newXpow_sum/newNoiseFloor);
%      figure(10);
%      subplot(1,3,1);
%      imagesc(Xpow_max_snr);set(gca,'YDir','normal');title("取最大");
%      subplot(1,3,2);
%      imagesc(Xpow_sum_snr);set(gca,'YDir','normal');title("取和");
%      subplot(1,3,3);
%      imagesc(newXpow_sum_snr);set(gca,'YDir','normal');title("取部分和");
%      Xpow = max(Xpow,[],3); % 这里跟我提供的代码不一样，他里面是sum起来，但是取max似乎是差不多的意义
     % 取一部分range-angle heat map
     Xpow1 = generateHiSNRRAmap(dopplerOut,Xpow);  % 注意，如果转换成log，差别变小太多，我建议不要转换直接算，会带来计算压力倒是，那其实原因就是FFT的信噪比不如MUSIC啥的高
     Xsnr = Xpow1;
     if Xpow1(1,1) ~= 0
         Xsnr = 20*log(Xpow1);
     end
     [rr,dd,aa] = size(Xpow);
     Xpow_mult = ones([rr,dd]);
     for ang = 1:aa
         Xpow_mult = Xpow_mult.*(20*log10(Xpow(:,:,ang)));
     end
%      Xpow_mult = 20*log10(Xpow_mult);
     figure(2);subplot(2,2,2);imagesc(agl_grid,distanceCoor,Xpow_mult);
     set(gca,'YDIR','normal');
     Xpow = 20*log10(sum(Xpow,3));
%      if Xpow > 0
%         Xsnr = 20*log10(Xpow);
%      end
%      Xsnr = pow2db(Xpow);
%      Xpow = squeeze(sum(Xpow,3)/size(Xpow,3));
%      Xsnr = pow2db(Xpow/noisefloor);
%      save Xsnr Xsnr
     XData = agl_grid;
     YData = distanceCoor;
%      CData = fliplr(Xsnr);
	 CData = Xsnr;
     % 转换为笛卡尔坐标系--杨炎龙
     % 毫米波雷达
%      [carteXSNR,x,y] = polar2carte(Xsnr, agl_grid, distanceCoor, 2);
%      figure;
%      subplot(1,2,1);
%      imagesc(x,y,fliplr(carteXSNR));
%      title("毫米波");
%      set(gca,'YDir','normal'); 
%      % 激光雷达读取数据
%      lidarFrameCnt = 326;  % 单5线数据
%      lidarFrameCnt = 109;  % 一帧帧数据
%      lidarData = load('data/data_20210414/lidar/05_横着走1.txt');
     frame = frameIndex;
%      lidarAngleGrid = (lidarData(lidarFrameCnt*(frame-1)+1:lidarFrameCnt*frame,1)*256+lidarData(lidarFrameCnt*(frame-1)+1:lidarFrameCnt*frame,2))/100 - 90;
%      lidarAngleGrid = (lidarData(lidarFrameCnt*frame-(frame-2):lidarFrameCnt*(frame+1)-frame,1)*256+lidarData(lidarFrameCnt*frame-(frame-2):lidarFrameCnt*(frame+1)-frame,2))/100 - 90;
%      lidarData = lidarData(lidarFrameCnt*frame-(frame-2):lidarFrameCnt*(frame+1)-frame,13:end); % 去掉每包前13个数
     lidarData = lidarDataFrame(1:end-1,:,frame);
     lidarAngleGrid = (lidarData(:,1)*256+lidarData(:,2))/100-90;
     startNum = 40;
     endNum = 270;
     lidarData = lidarData(:,startNum:end);  % 去掉每包前13个数
     [m,n] = size(lidarData);
     lidarRangeGrid = 0.15*([1:n]);
%      [carteLidarXSNR,x_lidar,y_lidar] = polar2carte(lidarData', lidarAngleGrid, lidarRangeGrid, 8);
%      subplot(1,2,2);
%      imagesc(x_lidar,y_lidar,carteLidarXSNR);
%      title("激光雷达");
%      set(gca,'YDir','normal'); 
%      2d cfar 提取目标点
%      biMatRAMap = CFAR2d(Xpow, 1.005);
%      fusionLidarADmap = fusion(Xsnr,agl_grid,distanceCoor,lidarData,lidarAngleGrid);
     [radarConfMap,biMatRAMap,radarMiuSig] = generateRadarConfMap(Xsnr,agl_grid,distanceCoor);
%      [lidarConfMap,bitMatLidarMap,lidarMiuSig] = generateLidarConfMap(lidarData,lidarAngleGrid,lidarRangeGrid);
     [lidarConfMap,lidarMiuSig] = generateLidarConfMap_v2(lidarData,lidarAngleGrid,lidarRangeGrid);
%      [outconfMap,radarInLidarConfMap,new_lidarAngleGrid,new_lidarRangeGrid] = fusionV0_1(radarConfMap,lidarConfMap,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,2);
     [confMap,allMiuSig,twoDimg,new_AngleGrid,new_RangeGrid] = fusionV0_2(radarMiuSig,lidarMiuSig,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid);
%      [carteConfMap,x_cate,y_cate] = polar2carte(confMap, new_AngleGrid, new_RangeGrid, 5);
%      [outFusionConfMap,radarRawInLidarMap,conf_angleGrid,conf_rangeGrid] = fusionV0_1(radarConfMap,lidarConfMap,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
%      radarMiuSig=size(radarMiuSig)
%      lidarMiuSig = size(lidarMiuSig)
     %% 直接把两张原始图片插一起
%      radarRawMap = pow2db(Xsnr);
     radarRawMap = Xsnr;
     maxRadar = max(radarRawMap(:));
     radarRawMap_norm = radarRawMap/maxRadar;
     if maxRadar <= 0
        radarRawMap_norm = zeros(size(radarRawMap_norm));
     end
     maxLidar = max(lidarData(:));
     lidarData_norm = lidarData/(maxLidar);
%      [outRawFusionMap,radarRawInLidarMap,new_angleGrid,new_rangeGrid] = fusionV0_1(radarRawMap_norm,lidarData_norm,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
     [outRawFusionMap,radarRawInLidarMap,new_angleGrid,new_rangeGrid] = fusionV0_1(Xpow_mult,lidarData_norm,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
     figure(4);
     subplot(1,2,1);imagesc(lidarAngleGrid,lidarRangeGrid,radarRawInLidarMap);
     set(gca,'YDIR','normal');
     subplot(1,2,2);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData_norm');
     set(gca,'YDIR','normal');
%      [straightConfMap,straightMiuSig] = generateLidarConfMap_v2(outRawFusionMap',new_angleGrid,new_rangeGrid);
     % 存一帧数据分析下特征
%      save outRawFusionMap.mat outRawFusionMap
%      figure(3);
%      colormap(jet);
%      subplot(2,3,1);imagesc(agl_grid,distanceCoor,radarConfMap);
%      set(gca,'YDIR','normal');title("毫米波置信度图");
%      subplot(2,3,2);imagesc(lidarAngleGrid,lidarRangeGrid,lidarConfMap');
%      set(gca,'YDIR','normal');title("激光置信度图");
% %      subplot(2,3,3);imagesc(new_angleGrid,new_rangeGrid,straightConfMap');
% %      set(gca,'YDIR','normal');title("叠加输入的置信度图");
%      subplot(2,3,4);imagesc(conf_angleGrid,conf_rangeGrid,outFusionConfMap);
%      set(gca,'YDIR','normal');title("叠加置信度图");
%      figure(4);subplot(2,2,1);imagesc(agl_grid,distanceCoor,Xsnr);
%      set(gca,'YDIR','normal');title("毫米波原始信号");
%      figure(2);
%      colormap(jet);
%      subplot(4,2,2);imagesc(new_AngleGrid,new_RangeGrid,twoDimg);
%      set(gca,'YDIR','normal');
%      subplot(4,2,[2,4]);imagesc(new_AngleGrid,new_RangeGrid,confMap);
%      set(gca,'YDIR','normal');title("融合置信度图");

%      radarInLidarConfMap = [];
%     endNum = 300;
%      subplot(4,2,[1,3]);imagesc(new_angleGrid,new_rangeGrid(1:end-endNum),outRawFusionMap(1:end-endNum,:));
%      subplot(4,2,[1,3]);imagesc(new_angleGrid,new_rangeGrid,outRawFusionMap);
%      a = outRawFusionMap(1:end-200,:);
%      subplot(4,2,[1,3]);imagesc(a);
%      set(gca,'YDIR','normal');title("输入");
%      %% 将静态毫米波图与激光雷达图插在一起
%      maxXpowRaw = max(XpowRaw(:));
%      XpowRaw_norm = XpowRaw/maxXpowRaw;
%      if maxXpowRaw <= 0
%          XpowRaw_norm = zeros(size(XpowRaw));
%      end
%      [outXpowRawFusionMap,XpowRawInLidarMap,~] = fusionV0_1(XpowRaw_norm,lidarData_norm,agl_grid,distanceCoor,lidarAngleGrid,lidarRangeGrid,1);
%      % 画合成图
%     dopplerLog2Abs = log2(abs(dopplerOut));
%     dopplerSum = sum(dopplerLog2Abs, [3 4]);
%     dopplerSum = squeeze(dopplerSum); 
%      figure(4);
% %     subplot(2,4,5);imagesc(velocityCoor,distanceCoor,dopplerSum);title('毫米波DopplerFFT波形');
% %     set(gca,'YDir','normal');
% %     xlabel('速度(m/s)');ylabel('距离(m)');zlabel('频谱幅值');
% %      subplot(2,4,5);imagesc(new_angleGrid,new_rangeGrid,XpowRawInLidarMap);
%      subplot(2,4,5);imagesc(XData,YData,20*log10(XpowRaw));
%      set(gca,'YDir','normal');
%      xlabel('角度');ylabel('距离(m)');title('毫米波带静态原始信号(20*log10(pow))');
%     subplot(2,4,6);imagesc(new_lidarAngleGrid,new_lidarRangeGrid,radarInLidarConfMap);
%     set(gca,'YDir','normal');
%     title("毫米波雷达confMap");
%     xlabel("角度");
%     ylabel("距离");
%     subplot(2,4,7);imagesc(lidarAngleGrid,lidarRangeGrid,lidarConfMap');
%     set(gca,'YDir','normal');
%     title("激光雷达confMap");
%     xlabel("角度");
%     subplot(2,4,8);imagesc(lidarAngleGrid,lidarRangeGrid,outconfMap);
%     set(gca,'YDir','normal');
%     title("融合");
%     xlabel("角度");
%      % 结束
     figure(2);subplot(2,2,3);imagesc(XData,YData,CData); % 
     set(gca,'YDir','normal');title("速度过滤");
     subplot(2,2,4);imagesc(XData,YData,Xpow); % 
     set(gca,'YDir','normal');title("不过滤");
%      title('毫米波动态原始信号');
% %      subplot(1,4,2);imagesc(XData(27:38),YData,fliplr(biMatRAMap(:,27:38))); % (:,27:38)
%      subplot(2,4,2);imagesc(XData,YData,biMatRAMap); % (:,27:38)
%      set(gca,'YDir','normal');
%      title('毫米波CFAR结果');
%      subplot(2,4,3);imagesc(lidarAngleGrid,lidarRangeGrid,lidarData');
%      set(gca,'YDir','normal');
%      title('激光雷达原始数据');
% %      subplot(1,4,4);plot(lidarData(54,1:224)');
%      subplot(2,4,4);imagesc(lidarAngleGrid,lidarRangeGrid,bitMatLidarMap');
%      set(gca,'YDir','normal');
%      title('激光雷达CFAR结果');
     
     %% 画笛卡尔坐标图
%     [carteRawFusionMap,x_cate,y_cate] = polar2carte(outRawFusionMap, new_angleGrid, new_rangeGrid, 20);
%     [carteRawFusionMap,x_cate,y_cate] = polar2carte(Xpow,agl_grid-5,distanceCoor, 8);
%     figure(2);
% %     clf()
%     subplot(2,2,2);
%     imagesc(x_cate,y_cate,carteRawFusionMap);
%      set(gca,'YDir','normal');
%      xlabel('x(m)');ylabel('y(m)');title('原始信号直接合成');
%      subplot(2,4,[3,7]);
%      imagesc(new_angleGrid,new_rangeGrid,outXpowRawFusionMap);
%      set(gca,'YDir','normal');
%      xlabel('角度');ylabel('距离(m)');title('激光融合毫米波(带静态)');
%      subplot(2,4,8);imagesc(new_AngleGrid,new_RangeGrid,confMap);
%      set(gca,'YDIR','normal');
% %      xlabel('角度');ylabel('距离(m)');title('新的置信度图');
    %%  存点mat
%     saveRAmaps(lidarData',radarRawInLidarMap,frameIndex,"./data/saveRAmap_0513/");
%     path = "./data/data_save/save7mfogfusion_0514/";
%     fileCnt = num2str(frameIndex);
%     fusion_filename = strcat(path,"fusion_",fileCnt,".mat");
%     save(fusion_filename,'outRawFusionMap');
    %% 显示点云
     [rows,cols] = size(allMiuSig);
     if rows~=0 || cols~=0
         % 由于线序是3，4，5
         if lineId <= 3
             lineId = lineId+2;
%              pcStrc = calElevAngLidar(lineId,allMiuSig);
         else
             lineId = lineId-3+2+8;
         end
        pcStrc = calElevAngLidar(lineId,allMiuSig);
%         angle = allMiuSig(:,1)*pi/180;
%         range = allMiuSig(:,2);
%         pcStrc.vertex.x = range.*sin(angle);
%         pcStrc.vertex.y = range.*cos(angle);
%         pcStrc.vertex.z = zeros(length(angle),1);
%         pcStrc.vertex.i = allMiuSig(:,3);
%         figure(3);
%         scatter3(pcStrc.vertex.x,pcStrc.vertex.y,pcStrc.vertex.z,50,'.');view([0,0,1]);
%         xlim([-10,10]);ylim([0,80]);title("点云");
%         grid off;axis on;
%         filename = sprintf('./data/data_20210607深大北门/mrr+urr-urr-ply/mrr+urr-urr%06d.ply',frameIndex);
%         pcwrite(pc, filename, 'Encoding', 'ascii');
%         PLY_WRITE( pcStrc, filename, 'ascii' );
     end
end