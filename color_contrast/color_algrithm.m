clear all; close all;
% 1.启发式染色算法将时隙划分成段
% 2.优先分配两跳以内被占用的时隙
% 3.分配时取两节点间最小的段大小，可将剩余部分扩大
% 评价标准
% 吞吐量
% 公平性

% 参数设置
% N = 16; % 节点数量
M = 50; % 时隙数量
Ne = 10; % 试验次数
index_N = 1;

for N = 5:35

throughput_all = 100;
throughput_var = 0;
for ll = 1 : Ne
clearvars -except N M Ne index_N throughput_all throughput_var throughput_all_rec throughput_var_rec 
table = zeros(N, M); % 时隙分配列表

% 生成网络拓扑
while true
    topo_tmp = randi([0,1], [N,N]);
    topo = topo_tmp + topo_tmp';
    topo(topo ~= 1) = 0;
    judge_tmp = sum(topo, 2);
    if isempty(find(judge_tmp == 0, 1))
        break;
    end
end

% [~,~,topo,~] = mknet(N,4,2,2);
% topo = double(topo) - eye(N);
% topo = [0 1 1 1 0; 1 0 0 0 0; 1 0 0 0 0; 1 0 0 0 1; 0 0 0 1 0];
% topo = [0 1 0 0; 1 0 1 0; 0 1 0 1; 0 0 0 1];

% topo = ones(N) - eye(N);  % 全连接网络

% 染色算法
color_table = zeros(N,2); % 每个节点邻接的颜色（未分配），剩余时隙数目
color_table(:,2) = M; % 初始化剩余时隙数目
% 计算邻接的颜色
% 找到邻居节点和邻居间的连接关系

link = cell(N,1);

for source = 1:N
    neigbor = [];
    for i = 1:N
        if topo(source, i)
            neigbor = [neigbor, i];
            link_tmp = link{source,1};
            link_tmp = [link_tmp;[source, i, 0]];
            link{source,1} = link_tmp;
        end
    end

    for i = 1:length(neigbor)
        for j = i+1:length(neigbor)
            if topo(neigbor(i), neigbor(j))
                link_tmp = link{source,1};
                link_tmp = [link_tmp;[neigbor(i), neigbor(j), 0]];
                link{source,1} = link_tmp;
            end
        end
    end

    link_color = link{source,1}(:,3);
    while ~isempty(find(link_color == 0, 1))
        color_cnt_max = -1;
        color_set_max = [];
        index = 0;
        link_tmp = link{source,1};
        for i = 1:size(link_tmp, 1)
            if link_tmp(i, 3) == 0
                [color_set, color_cnt] = cnt_color(link_tmp(i, 1), link_tmp(i, 2), link_tmp);
                if color_cnt > color_cnt_max
                    color_cnt_max = color_cnt;
                    color_set_max = color_set;
                    index = i;
                end
            end
        end
        if isempty(color_set_max)
            link_tmp(index, 3) = 1;
            link{source,1} = link_tmp;
            link_color = link{source,1}(:,3);
            continue;
        end
        flag = 0;
        for l = 1:max(color_set_max)
            if isempty(find(color_set_max == l, 1))
                link_tmp(index, 3) = l;
                flag = 1;
                break;
            end
        end
        if flag == 0
            link_tmp(index, 3) = l+1;
        end
        link{source,1} = link_tmp;
        link_color = link{source,1}(:,3);
    end
    color_table(source, 1) = max(link_color);
end

% 进行时隙分配
% 随机挑选链路进行
k = 1;
for i = 1:N
    for j = i+1:N
        if topo(i,j)
            link_t{k} = [i,j];
            k = k + 1;
        end
    end
end

mark = zeros(1, length(link_t));
num_link_all = zeros(size(mark));
num_link_index = 1;

while sum(mark) < length(link_t)
    index_link = find(mark == 0);
    sel = randi([1, length(index_link)]);
    mark(index_link(sel)) = 1;
    [table, color_table, num_link] = allocate(link_t{index_link(sel)}(1), link_t{index_link(sel)}(2), table, color_table, topo);
    num_link_all(num_link_index) = num_link;
    num_link_index = num_link_index + 1;
end

% 性能评估
throughput = zeros(N,1);
for i = 1:N
    throughput(i) = length(find(table(i,:) ~= 0));
end

throughput_all = min(throughput_all, sum(throughput)/2/M);
throughput_var = throughput_var + var(num_link_all)/M;

end

throughput_all_rec(index_N) = throughput_all;
throughput_var_rec(index_N) = throughput_var/Ne;
index_N = index_N + 1;

end

figure
plot(5:35, throughput_all_rec)

figure
plot(5:35, throughput_var_rec)

[1.34000000000000,1.84000000000000,2.02000000000000,2.54000000000000,2.52000000000000,3.30000000000000,3.68000000000000,4.04000000000000,4.72000000000000,5.28000000000000,5.38000000000000,6.32000000000000,6.62000000000000,7.10000000000000,7.54000000000000,7.84000000000000,8.42000000000000,8.86000000000000,9.38000000000000,9.98000000000000,10.3000000000000,10.8200000000000,11.4400000000000,12.0400000000000,12.0400000000000,12.8800000000000,13.3400000000000,13.8200000000000,14.4200000000000,14.8400000000000,15.3400000000000]
