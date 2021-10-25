%% 连续读取数据

    fusionConfMap = fogdetector_v1("./data/data_save/saveYuanliSZ_0514_net/",5);
    fusionConfMap = fogdetector_v1("./data/data_save/save7mfogfusion_0514/",5);
    figure(10);imagesc(lidarRAmap);set(gca,'YDir','normal');
    a = dir("./data/saveFusionMap_0514/");
    a.name
    