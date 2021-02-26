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

for N = 16

throughput_all = 0;
throughput_var = 0;
for ll = 1 : Ne

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

% topo = ones(N) - eye(N);  % 全连接网络

% 染色算法
color_table = zeros(N,2); % 每个节点邻接的颜色（未分配），剩余时隙数目
color_table(:,2) = M; % 初始化剩余时隙数目
color_table(:,1) = sum(topo, 2);

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

throughput_all = throughput_all + sum(throughput)/2/M;
throughput_var = throughput_var + var(num_link_all)/M;

end

throughput_all_rec(index_N) = throughput_all/Ne;
throughput_var_rec(index_N) = throughput_var/Ne;
index_N = index_N + 1;

end

% figure
% plot(5:35, throughput_all_rec)

% figure
% plot(5:35, throughput_var_rec)

% 1 hop 邻居最小数目
% N/(min(sum(topo, 2))+1)

