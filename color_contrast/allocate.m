function [table, color_table, num_allocate_fnl] = allocate(src, dst, table, color_table, topo)
%allocate - Description
%
% Syntax: [table] = allocate(src, dst, table, color_table, topo)
%
% Long description
num_allocate = min(ceil(color_table(src, 2) / color_table(src, 1)), ceil(color_table(dst, 2) / color_table(dst, 1)));
% 计算两跳邻居
two_hop_src = [];
two_hop_dst = [];
N = size(topo, 1);
M = size(table, 2);
for i = 1:N
    if topo(src, i)
        for j = 1:N
            if topo(i, j) && topo(src, j) == 0 && length(find(two_hop_src == j)) == 0
                two_hop_src = [two_hop_src, j];
            end
        end
    end
end

for i = 1:N
    if topo(dst, i)
        for j = 1:N
            if topo(i, j) && topo(dst, j) == 0 && length(find(two_hop_dst == j)) == 0
                two_hop_dst = [two_hop_dst, j];
            end
        end
    end
end

% 统计可分配的时隙
slot_tmp = table(src, :) + table(dst, :);
slot_free = find(slot_tmp == 0);
slot_two_hop_src_occupy = [];
slot_two_hop_dst_occupy = [];

for i = 1:length(two_hop_src)
    slot_two_hop_src_occupy = [slot_two_hop_src_occupy, find(table(two_hop_src(i), :) ~= 0)];
end

for i = 1:length(two_hop_dst)
    slot_two_hop_dst_occupy = [slot_two_hop_dst_occupy, find(table(two_hop_dst(i), :) ~= 0)];
end

slot_two_hop_src_occupy = unique(slot_two_hop_src_occupy);
slot_two_hop_dst_occupy = unique(slot_two_hop_dst_occupy);

% for i = 1:length(slot_free)
%     if length(find(slot_two_hop_src_occupy == mod(slot_free(i)-1, M)+1))
%         slot_free(i) = slot_free(i) - M;
%     end
%     if length(find(slot_two_hop_dst_occupy == mod(slot_free(i)-1, M)+1))
%         slot_free(i) = slot_free(i) - M;
%     end
% end

% slot_free = sort(slot_free);
% 按优先级分配时隙
num_allocate_fnl = min(length(slot_free), num_allocate);
rand_index = randperm(length(slot_free));
slot_free = slot_free(rand_index);
for i = 1:num_allocate_fnl
    index = mod(slot_free(i)-1, M)+1;
    table(src, index) = dst;
    table(dst, index) = src;
end
    
% 修改color_table
% color_table(src, 1) = color_table(src, 1) - 1;
% color_table(src, 2) = color_table(src, 2) - num_allocate_fnl;

% color_table(dst, 1) = color_table(dst, 1) - 1;
% color_table(dst, 2) = color_table(dst, 2) - num_allocate_fnl;

end