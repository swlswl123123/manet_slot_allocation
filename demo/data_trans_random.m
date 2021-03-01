function [random_trans, random_delay] = data_trans_random(N, M, R, t_last, topo)
% 1. 按链路与时间均匀生成数据包
% 2. 计算转发表，原始列表发送/接收交替
% 3. 每个节点维护传输邻居节点个数的传输队列cell([src, dst, pkid])，
% 4. 每一帧开始时按负载量决定实际使用的时隙数目，发送方申请（轮流抢占，其余不超过原始列表），随机链路，轮流两节点优先级
% 5. 按时间片遍历，传送负载
% 参数设置
% N = 5; % 节点数量
% M = 50; % 时隙数量

table = zeros(N, M); % 时隙分配列表

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               完成预分配
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 生成网络拓扑
% while true
%     topo_tmp = randi([0,1], [N,N]);
%     topo = topo_tmp + topo_tmp';
%     topo(topo ~= 1) = 0;
%     judge_tmp = sum(topo, 2);
%     if isempty(find(judge_tmp == 0, 1))
%         break;
%     end
% end

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
    table = allocate_random(link_t{index_link(sel)}(1), link_t{index_link(sel)}(2), table, 1);
    table = allocate_random(link_t{index_link(sel)}(2), link_t{index_link(sel)}(1), table, 1);
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

% t_last = 10;
% R = 1000;
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

% 按时间扫描
% 每个时隙为1ms，一帧为50ms
cur_time = 0;
flag = 1;
while cur_time < t_last*2
    % 每帧开始时刻 
    % 计算实际需求的时隙
    slot_need_table = zeros(N, N);
    table = table_origin;
    for i = 1:N
        for j = 1:N
            if topo(i, j)
                slot_need_table(i, j) = ceil(size(txQ{i, j}, 1) / pk_rate);
            end
        end
    end
    % 分配扩展时隙
    mark = zeros(1, length(link_t));
    while sum(mark) < length(link_t)
        index_link = find(mark == 0);
        sel = randi([1, length(index_link)]);
        mark(index_link(sel)) = 1;
        ext_slot_src = link_t{index_link(sel)}(1);
        ext_slot_dst = link_t{index_link(sel)}(2);
        if flag % src主导
            num = slot_need_table(ext_slot_src, ext_slot_dst);
            table = allocate_random(ext_slot_src, ext_slot_dst, table, num);
        else    % dst主导
            num = slot_need_table(ext_slot_dst, ext_slot_src);
            table = allocate_random(ext_slot_dst, ext_slot_src, table, num);
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
            if table(j, i) > 0 && table(table(j, i), i) == -j % 收发对准
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
end

delay = recv_t - send_t;

% figure
% plot(time(1:t_last/0.01), throughput(2:t_last/0.01+1))
random_trans = throughput(2:t_last/0.01+1);
random_delay = delay(delay > 0);

% save('random_trans.mat', 'random_trans');
% save('random_delay.mat', 'random_delay');

% figure
% plot(pk_send)

% drop 大于 10 % 便不统计
end