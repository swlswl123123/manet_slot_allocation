clear all; close all;
% rng(0)
N = 16;
M = 50;
R = 1000;
t_last = 10;

while true
    topo_tmp = randi([0,1], [N,N]);
    topo = topo_tmp + topo_tmp';
    topo(topo ~= 1) = 0;
    judge_tmp = sum(topo, 2);
    if isempty(find(judge_tmp == 0, 1))
        break;
    end
end

[advanced_trans, advanced_delay] = data_trans_advanced_fix(N, M, R, t_last, topo);
[color_trans, color_delay] = data_trans_color_fix(N, M, R, t_last, topo);
[random_trans, random_delay] = data_trans_random(N, M, R, t_last, topo);

figure
plot(advanced_trans, '-o')
hold on;
plot(color_trans, '-*')
hold on;
plot(random_trans, '-^')

legend('advanced', 'colored', 'random')

mean(advanced_delay)
mean(color_delay)
mean(random_delay)
