function [match_right] = KMhungarian(adj_matrix, N)

global N
global adj_matrix
global label_left
global label_right
global match_right
global visit_left
global visit_right
% 左右各有N个点
% KM算法要求左右两边的节点数相等，可以通过添加虚拟节点的方法实现
% N = 5;
% adj_matrix = [3 4 6 4 9;
%     6 4 5 3 8;
%     7 5 3 4 2;
%     6 3 2 2 5;
%     8 4 5 4 7];
% N = 9;
% adj_matrix = round(rand(N)*100);
% adj_matrix(adj_matrix<70) = 0;
% 初始化顶标
label_left = max(adj_matrix, [], 2);
label_right = zeros(N, 1);
% 初始化匹配结果
match_right = ones(N, 1) * nan;
% 初始化辅助变量
visit_left = ones(N, 1) * false;
visit_right = ones(N, 1) * false;
res = KM();
end

% KM主函数
function res = KM()
global N
global adj_matrix
global label_left
global label_right
global match_right
global visit_left
global visit_right

graph_num = 1;
% display_graph(graph_num, '原始二部图');
% 对左边的点依次进行处理
for i = 1: N
    while 1
        % 重置辅助变量
        visit_left = ones(N, 1) * false;
        visit_right = ones(N, 1) * false;
        % 能找到可行匹配
        if find_path(i)
            break;
        end
        % 不能找到可行匹配，修改顶标
        % (1)将所有在增广路中的X方点的label全部减去一个常数d
        % (2)将所有在增广路中的Y方点的label全部加上一个常数d
        d = Inf;
        for j = 1: N
            if visit_left(j)
               for k = 1: N
                   if ~visit_right(k)
                       % 左边的点中已经访问过的点，即已经匹配过的点可能需要重新匹配以得到更大的总权值，
                       % 所以修改顶标，往子图中添加一条边，重新寻找增广路看能不能增广
                       % 取与左边的点相邻的未匹配边中跟当前存在子图中的以该点为端点的边相差最小的两条边
                       % 这样才能保持总权值最大
                       d = min(d, label_left(j) + label_right(k) - adj_matrix(j, k));
                   end
               end
            end
        end
        for k = 1: N
            if visit_left(k)
                label_left(k) = label_left(k) - d;
            end
            if visit_right(k)
                label_right(k) = label_right(k) + d;
            end
        end
    end
    graph_num = graph_num + 1;
%     display_graph(graph_num, ['第' num2str(i) '步']);
end

graph_num = graph_num + 1;
display_graph(graph_num, '最终结果');

res = 0;
for j = 1: N
    if match_right(j) >=0 && match_right(j) < N
        res = res + adj_matrix(match_right(j), j);
    end
end
end

% 寻找增广路，深度优先
function result = find_path(i)
global N
global adj_matrix
global label_left
global label_right
global match_right
global visit_left
global visit_right
visit_left(i) = true;
for j = 1: length(adj_matrix(i, :))
    match_weight = adj_matrix(i, j);
    if visit_right(j)
        % 已被匹配（解决递归中的冲突）
        continue;
    end
    gap = label_left(i) + label_right(j) - match_weight;
    % 当 gap == 0 时 x_i 和 y_j 之间存在一条边，且该边是当前 x_i 可以匹配的权值最大的边
    if gap == 0
        % 找到可行匹配
        visit_right(j) = true;
        % j未被匹配，或虽然j已被匹配，但是j的已匹配对象有其他可选备胎
        % 此处同匈牙利算法
        if isnan(match_right(j)) || find_path(match_right(j))
            match_right(j) = i;
            result = true;
            return;
        end
    end
end
result = false;
return;
end


function display_graph(graph_num, title_name)
global N
global adj_matrix
global label_left
global label_right
global match_right
global visit_left
global visit_right
figure(graph_num);
set(gcf, 'Position', [100, 100, 1000, 500]);
subplot(1,2,1);
cla
set( gca, 'XTick', [], 'YTick', [] );
set( gca, 'TickLength', [0 0]);
box on
hold on
xlim([-1, 2]);
ylim([-N, 1]);
temp = N - 1;
% 绘制匹配边
for j = 1: length(match_right)
    if ~isnan(match_right(j))
        plot([0, 1], [(match_right(j) - 1) * -temp/(N-1), (j - 1) * -temp/(N-1)], 'r', 'LineWidth', 2); 
    end
end
% 绘制点
scatter(zeros(1, N), 0:-temp/(N-1):-temp, 20, [217/255 83/255 25/255], 'filled');
scatter(ones(1, N), 0:-temp/(N-1):-temp, 20, [0 114/255 189/255], 'filled');
for j = 1: N
    text(-0.2, (j - 1) * -temp/(N-1), ['x_' num2str(j)]);
end
for j = 1: N
    text(1.1, (j - 1) * -temp/(N-1), ['y_' num2str(j)]);
end
for j = 1: N
    text(-0.5, (j - 1) * -temp/(N-1), num2str(label_left(j)), 'Color', [217/255 83/255 25/255], 'FontSize', 15);
end
for j = 1: N
    text(1.4, (j - 1) * -temp/(N-1), num2str(label_right(j)), 'Color', [0 114/255 189/255], 'FontSize', 15);
end

title(title_name);

% 边权重
subplot(1,2,2);
cla
set( gca, 'XTick', [], 'YTick', [] );
set( gca, 'TickLength', [0 0]);
box on
hold on
xlim([-1, N]);
ylim([-N, 1]);
for j = 0 : N-1
   plot([j, j], [-N, 1], 'k') ;
end
for j = -N+1 : 0
   plot([-1, N], [j, j], 'k') ;
end
for j = 1: N
    display_text(j, 0, ['x_' num2str(j)], 0);
end
for j = 1: N
    display_text(0, j, ['y_' num2str(j)], 0);
end
for i = 1: N
    for j = 1: N
        if match_right(j) == i
            display_text(i, j, num2str(adj_matrix(i,j)), 1);
        else
            display_text(i, j, num2str(adj_matrix(i,j)), 0);
        end
    end
end
% saveas(gcf, [num2str(graph_num) '.png']);
end

% 在x行y列显示t文本
function display_text(x, y, t, bold)
if bold
    text(y - 0.5, -x + 0.5, t, 'FontSize', 15, 'FontWeight', 'bold', 'Color', 'r');
else
    text(y - 0.5, -x + 0.5, t, 'FontSize', 15, 'FontWeight', 'normal', 'Color', 'k');
end
end