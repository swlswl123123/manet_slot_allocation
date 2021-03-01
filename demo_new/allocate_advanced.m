function [table, color_table, num_allocate_fnl] = allocate_advanced(src, dst, table, color_table, topo, type)
    %allocate - Description
    %
    % Syntax: [table] = allocate(src, dst, table, color_table, topo)
    %
    % Long description
    N = size(topo, 1);
    M = size(table, 2);
    if type == -1    % 分配基本时隙
        slot_tmp = abs(table(src, :)) + abs(table(dst, :));
        slot_free = find(slot_tmp == 0);
        num_allocate_fnl = min(length(slot_free), 1);
        % 随机打乱
        slot_free = slot_free(randperm(length(slot_free)));
        for i = 1:num_allocate_fnl
            index = mod(slot_free(i)-1, M)+1;
            table(src, index) = dst;
            table(dst, index) = -src;
        end
        color_table(src, 2) = color_table(src, 2) - num_allocate_fnl;  
        color_table(dst, 2) = color_table(dst, 2) - num_allocate_fnl; 
    else
        num_allocate = type;
        % 计算两跳邻居
        two_hop_src = [];
        two_hop_dst = [];
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
        slot_tmp = abs(table(src, :)) + abs(table(dst, :));
        slot_free = find(slot_tmp == 0);

        % 统计两跳邻居已占用的时隙
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
        
        % 对于被占用的时隙赋予较高的优先级
        for i = 1:length(slot_free)
            if length(find(slot_two_hop_src_occupy == mod(slot_free(i)-1, M)+1))
                slot_free(i) = slot_free(i) - M;
            end
            if length(find(slot_two_hop_dst_occupy == mod(slot_free(i)-1, M)+1))
                slot_free(i) = slot_free(i) - M;
            end
        end
        
        % 按优先级分配时隙
        slot_free = sort(slot_free);
        num_allocate_fnl = min(length(slot_free), num_allocate);
        % 可考虑打乱
        prio_idx = length(find(slot_free < 0));
        slot_free_sub = slot_free(prio_idx+1:end);
        % 打乱顺序
        slot_free_sub = slot_free_sub(randperm(length(slot_free_sub)));
        slot_free = [slot_free(1:prio_idx), slot_free_sub];
        for i = 1:num_allocate_fnl
            index = mod(slot_free(i)-1, M)+1;
            table(src, index) = dst;
            table(dst, index) = -src;
        end
            
        % 修改color_table
        color_table(src, 1) = color_table(src, 1) - 1;
        color_table(src, 2) = color_table(src, 2) - num_allocate_fnl;    
        color_table(dst, 1) = color_table(dst, 1) - 1;
        color_table(dst, 2) = color_table(dst, 2) - num_allocate_fnl;
    end
end