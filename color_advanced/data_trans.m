clear all; close all;
% 1. 按链路与时间均匀生成数据包
% 2. 计算转发表，原始列表发送/接收交替
% 3. 每个节点维护传输邻居节点个数的传输队列cell([src, dst, pkid])，
% 4. 每一帧开始时按负载量决定实际使用的时隙数目，发送方申请（轮流抢占，其余不超过原始列表），随机链路，轮流两节点优先级
% 5. 按时间片遍历，传送负载
rng(1)
% 参数设置
N = 16; % 节点数量
M = 50; % 时隙数量

table = zeros(N, M); % 时隙分配列表

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               完成预分配
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% 先分配基本时隙
mark = zeros(1, length(link_t));

while sum(mark) < length(link_t)
    index_link = find(mark == 0);
    sel = randi([1, length(index_link)]);
    mark(index_link(sel)) = 1;
    [table, color_table, num_link] = allocate(link_t{index_link(sel)}(1), link_t{index_link(sel)}(2), table, color_table, topo, 0);
end

% 分配扩展时隙
mark = zeros(1, length(link_t));
num_link_all = zeros(size(mark));
num_link_index = 1;

while sum(mark) < length(link_t)
    index_link = find(mark == 0);
    sel = randi([1, length(index_link)]);
    mark(index_link(sel)) = 1;
    [table, color_table, num_link] = allocate(link_t{index_link(sel)}(1), link_t{index_link(sel)}(2), table, color_table, topo, 1);
    num_link_all(num_link_index) = num_link+2;
    num_link_index = num_link_index + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               计算路由表
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forword_table = cell(16,1);
[p,q,v]=find(topo);
topo_s=sparse(p,q,v,N,N);
for src_idx = 1:N
    [~,path,~] = graphshortestpath(topo_s, src_idx);
    pred = zeros(1, N);
    for i = 1:N
        if i == src_idx
            continue;
        end
        pred(i) = path{i}(2);
    end
    forword_table{src_idx} = pred;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               生成负载
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table = table(:, randperm(M));
table_origin = table;

% 每包100 bits
% 物理层 1Mbps
% 10ms能发10包
pk_rate = 10;

t_last = 10;
R = 3000;
Np = R*N*t_last; % 包的数量

L = npermutek(1:N,2); % 所有链路
Li = randi(length(L), Np, 1);

id = (1:Np)';
src = L(Li,1);
dst = L(Li,2);
recv_t = zeros(Np, 1);
send_t = t_last*rand(Np,1);
time = [];
throughput = [0];
pk_send = [0];

% 传输队列
txQ = cell(N, N);
frame_cnt = 1; % 轮流计数

% 按时间扫描
% 每个时隙为1ms，一帧为50ms
cur_time = 0;
flag = 1;
while cur_time < t_last+4
    % 每帧开始时刻 
    % 计算实际需求的时隙
    if frame_cnt > N
        frame_cnt = 1;
    end
    table_real = zeros(N, M);
    table = table_origin;
    for i = 1:N
        slot_need_ext = 0;
        slot_ext_id = 0;
        for j = 1:N
            if topo(i, j)
                slot_need = ceil(size(txQ{i, j}, 1) / pk_rate);
                slot_pre = find(table(i, :) == j);
                slot_pre = slot_pre(randperm(length(slot_pre))); % 随机均匀
                if length(slot_pre) >= slot_need
                    slot_pre = slot_pre(1:slot_need);
                else
                    if slot_need - length(slot_pre) > slot_need_ext
                        slot_need_ext = slot_need - length(slot_pre);
                        slot_ext_id = j;
                    end
                end
                % slot_pre = slot_pre(1:min(slot_need, length(slot_pre)));
                table_real(i, slot_pre) = j;
            end
        end
        % 轮流抢占
        if frame_cnt == i && slot_need_ext ~= 0
            % 找到没被占用的部分
            ext_free = find(table_real(i, :) == 0);
            ext_free = ext_free(randperm(length(ext_free))); % 随机打乱
            ext_free = ext_free(1 : min(slot_need_ext, length(ext_free)));
            table_real(i, ext_free) = slot_ext_id;
            % for kk = 1:length(ext_free)
            %     table(table(i, ext_free(kk)), ext_free(kk)) = 0;
            %     table(table(slot_ext_id, ext_free(kk)), ext_free(kk)) = 0;
            % end
            table(i, ext_free) = slot_ext_id;
            table(slot_ext_id, ext_free) = i;
            if i < slot_ext_id
                flag = 1;
            else
                flag = 0;
            end
        end

    end
    % txQ_tmp = txQ;
    % table_real_tmp = table_real;
    % 链路交替优先分配时隙
    % 有问题：覆盖 收发匹配+-1不匹配
    for i = 1:N
        for j = i+1:N
            if topo(i, j)
                if flag % i占主导
                    % 统计j申请发送的时隙数并清零
                    j_req_send = length(find(table_real(j, :) == i));
                    table_real(j, find(table_real(j, :) == i)) = 0;
                    % 先满足i的申请
                    i_req_send = length(find(table_real(i, :) == j));
                    table_real(j, find(table_real(i, :) == j)) = -i;
                    % 再统计j中剩余的与i链路
                    ij_total = length(find(table(i, :) == j));
                    j_remain = find(table(j, :) == i & table_real(j, :) == 0);
                    j_remain = j_remain(randperm(length(j_remain)));
                    j_remain = j_remain(1 : min(j_req_send, ij_total - i_req_send));
                    % 满足j的申请
                    table_real(j, j_remain) = i;
                    table_real(i, j_remain) = -j;
                else % j占主导
                    % 统计i申请发送的时隙数并清零
                    i_req_send = length(find(table_real(i, :) == j));
                    table_real(i, find(table_real(i, :) == j)) = 0;
                    % 先满足j的申请
                    j_req_send = length(find(table_real(j, :) == i));
                    table_real(i, find(table_real(j, :) == i)) = -j;
                    % 再统计i中剩余的与j链路
                    ij_total = length(find(table(i, :) == j));
                    table_tmp = table_real(i, :) + table(i, :);
                    i_remain = find(table(i, :) == j & table_real(i, :) == 0);
                    i_remain = i_remain(randperm(length(i_remain)));
                    i_remain = i_remain(1 : min(i_req_send, ij_total - j_req_send));
                    % 满足i的申请
                    table_real(i, i_remain) = j;
                    table_real(j, i_remain) = -i;
                end
            end
        end
    end
    flag = 1 - flag;
    for i = 1:M
        % 放入队列中
        pk_idx = find(send_t >= cur_time & send_t < cur_time + 0.01);
        for p = 1:length(pk_idx)
            src_id = src(pk_idx(p));
            dst_id = dst(pk_idx(p));
            next_id = forword_table{src_id}(dst_id);
            q_tmp = txQ{src_id, next_id};
            q_tmp = [q_tmp;[src_id dst_id id(pk_idx(p))]];
            txQ{src_id, next_id} = q_tmp;
        end
        pk_send = [pk_send, pk_send(end) + length(pk_idx)];
        % 发送过程
        pk_trans = cell(N, 1);
        for j = 1:N
            if table(j, i) == 0
                continue;
            end
            if table_real(j, i) > 0 && table_real(table_real(j, i), i) == -j % 收发对准
                send_dst = table(j, i);
                q_tmp = txQ{j, send_dst};
                % 在本地删除
                if size(q_tmp, 1) > pk_rate
                    pk_trans{j} = q_tmp(1:pk_rate, :);
                    q_tmp = q_tmp(pk_rate+1:end, :);
                    txQ{j, send_dst} = q_tmp;
                else
                    txQ{j, send_dst} = [];
                    pk_trans{j} = q_tmp;
                end
            end
        end
        cur_time =  cur_time + 0.01;
        time = [time, cur_time];
        % 在目的地加入合适发送队列
        pk_send_num = 0;
        for j = 1:N
            if table(j, i) == 0
                continue;
            end
            trans_tmp = pk_trans{j};
            recv_dst = table(j, i);
            pk_send_num = pk_send_num + size(trans_tmp, 1);
            for k = 1:size(trans_tmp, 1)
                pk_dst = trans_tmp(k, 2);
                % 查路由表
                next_id = forword_table{recv_dst}(pk_dst);
                if forword_table{recv_dst}(pk_dst) == 0
                    recv_t(trans_tmp(k, 3)) = cur_time;
                else
                    q_tmp = txQ{recv_dst, next_id};
                    q_tmp = [q_tmp;trans_tmp(k, :)];
                    txQ{recv_dst, next_id} = q_tmp;
                end
            end
        end
        throughput = [throughput, (throughput(end)*(cur_time-0.01)+pk_send_num)/cur_time];
        % throughput = [throughput, throughput(end)+pk_send_num];
    end
    frame_cnt = frame_cnt + 1;
end

delay = recv_t - send_t;

figure
plot(time(1:t_last/0.01), throughput(2:t_last/0.01+1))
advanced_trans = throughput(2:t_last/0.01+1);
advanced_delay = delay(delay > 0);

save('advanced_trans.mat', 'advanced_trans');
save('advanced_delay.mat', 'advanced_delay');

% figure
% plot(pk_send)

% drop 大于 10 % 便不统计
