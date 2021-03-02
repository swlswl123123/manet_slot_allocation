clear all; close all;
% rng(0)
% 评价
% 超过某一时间认为loss，计算比较loss概率
% 如何找到极限吞吐量
% color为simplied的方法

throughput_advanced_rec = [];
throughput_color_rec = [];
throughput_random_rec = [];
delay_advanced_rec = [];
delay_color_rec = [];
delay_random_rec = [];

M = 50;
R = 100;
t_last = 5;
Ne = 1;

for N = 16
    throughput_advanced = 0;
    throughput_color = 0;
    throughput_random = 0;
    delay_advanced = 0;
    delay_color = 0;
    delay_random = 0;
    
    for ll = 1:Ne
        while true
            topo_tmp = randi([0,1], [N,N]);
            topo = topo_tmp + topo_tmp';
            topo(topo ~= 1) = 0;
            [p,q,v]=find(topo);
            topo_s=sparse(p,q,v,N,N);
            [dist,~,~] = graphshortestpath(topo_s, 1);
            if isempty(find(dist == inf, 1))
                break;
            end
        end
        
        [advanced_trans_fix, advanced_delay_fix] = data_trans_advanced_fix(N, M, R, t_last, topo);
        % [advanced_trans, advanced_delay] = data_trans_advanced(N, M, R, t_last, topo);
        [color_trans, color_delay] = data_trans_color_fix(N, M, R, t_last, topo);
        [random_trans, random_delay] = data_trans_random(N, M, R, t_last, topo);
        
        % figure
        % plot(advanced_trans_fix, '-o')
        % % hold on;
        % % plot(advanced_trans, '-d')
        % hold on;
        % plot(color_trans, '-*')
        % hold on;
        % plot(random_trans, '-^')
        
        % legend('advanced', 'colored', 'random')
        throughput_advanced = throughput_advanced + advanced_trans_fix(end);
        throughput_color = throughput_color + color_trans(end);
        throughput_random = throughput_random + random_trans(end);
        
        delay_advanced = delay_advanced + mean(advanced_delay_fix);
        delay_color = delay_color + mean(color_delay);
        delay_random = delay_random + mean(random_delay);
    end
    throughput_advanced_rec = [throughput_advanced_rec, throughput_advanced/Ne];
    throughput_color_rec = [throughput_color_rec, throughput_color/Ne];
    throughput_random_rec = [throughput_random_rec, throughput_random/Ne];
    delay_advanced_rec = [delay_advanced_rec, delay_advanced/Ne];
    delay_color_rec = [delay_color_rec, delay_color/Ne];
    delay_random_rec = [delay_random_rec, delay_random/Ne];
end

N_set = 5 : 16;
figure
plot(N_set, throughput_advanced_rec);
hold on;
plot(N_set, throughput_color_rec);
hold on;
plot(N_set, throughput_random_rec);

figure;
plot(N_set, delay_advanced_rec);
hold on;
plot(N_set, delay_color_rec);
hold on;
plot(N_set, delay_random_rec);


