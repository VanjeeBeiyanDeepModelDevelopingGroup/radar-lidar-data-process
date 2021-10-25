function [outImage] = project2image(targetMat,image_filepath)
    outImage = [];
    %% 将毫米波雷达或激光雷达检测到的物体投影到图像上去
    % 图像内参矩阵是多少
    k=[1155.023180,0.000000,967.178167;0.000000,1149.531704,541.645587;0.000000,0.000000,1.000000];
    % 外参
    % 毫米波到相机
    theta = 0/180*pi;
    Rr2c = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1]; % 旋转矩阵(绕z轴旋转)
    Tr2c = [0.27,-0.13,0]'; % 米
%     % 激光雷达到相机
%     theta = 0;
%     Rr2c = [cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1]; % 旋转矩阵(绕z轴旋转)
%     Tr2c = [270,0,0]; % 毫米
    % 读取图像
    image = imread(image_filepath);
    [rows,cols] = size(image(:,:,1));
    figure;imshow(image);
    %% 将confMap转换到笛卡尔坐标系来，再进行投影操作
    [numTarget,targetCols] = size(targetMat);
    uv_mat = [];
    xyz_mat = [];
    xyz_cam_mat = [];
    targetImage = [];
    targetMat
    for i=1:numTarget
        fprintf("=======================\n");
        range = targetMat(i,1);
        angle = targetMat(i,2)/180*pi;
        power = targetMat(i,3);
        x = range*sin(angle);
        y = 1.0;
        z = range*cos(angle);
        % 存所有的xyz
        xyz_mat(i,1) = x;
        xyz_mat(i,2) = y;
        xyz_mat(i,3) = z;
        xyz_mat(i,4) = power;
        fprintf("x,y,z: %f,%f,%f\n",x,y,z);
        % 移动 [x',y',z'] = R[x,y,z]+t
        pos_c = Rr2c*[x;y;z]+Tr2c;
        xyz_cam_mat(i,1) = pos_c(1);
        xyz_cam_mat(i,2) = pos_c(2);
        xyz_cam_mat(i,3) = pos_c(3);
        xyz_cam_mat(i,4) = power;
        % 投影
        uv_center = k*pos_c;
        uv_center = uv_center/uv_center(3);
        u = floor(uv_center(1)+cols/2);
        v = floor(rows/2-uv_center(2));
        if(u<0 || u > 10000 || v < 0 || v > 10000)
%         if(u<0 || v < 0 || u==inf || v == inf)
            continue
        end
        hold on;
        plot(u,v,'ro');
        text(u,v,num2str(power));
        % 存所有的uv
        uv_mat(i,1) = u;
        uv_mat(i,2) = v;
        uv_mat(i,3) = power;
        targetImage(u,v) = 1;
        fprintf("u,v: %d,%d\n",u,v);
    end
    figure;
    plot(xyz_mat(:,1),xyz_mat(:,3),'ro');
    text(xyz_mat(:,1),xyz_mat(:,3),num2str(xyz_mat(:,4)));
%     hold on
%     plot(xyz_cam_mat(:,1),xyz_cam_mat(:,3),'bo');
%     text(xyz_cam_mat(:,1),xyz_cam_mat(:,3),num2str(xyz_cam_mat(:,4)));
end

